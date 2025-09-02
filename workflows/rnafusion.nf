/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BUILD_REFERENCES              }   from '../subworkflows/local/build_references'
include { CAT_FASTQ                     }   from '../modules/nf-core/cat/fastq/main'
include { FASTQ_FASTQC_UMITOOLS_FASTP   }   from '../subworkflows/nf-core/fastq_fastqc_umitools_fastp/main'
include { QC_WORKFLOW                   }   from '../subworkflows/local/qc_workflow'
include { STARFUSION_DETECT             }   from '../modules/nf-core/starfusion/detect/main'
include { STRINGTIE_WORKFLOW            }   from '../subworkflows/local/stringtie_workflow/main'
include { FUSIONCATCHER_WORKFLOW        }   from '../subworkflows/local/fusioncatcher_workflow'
include { FUSIONINSPECTOR_WORKFLOW      }   from '../subworkflows/local/fusioninspector_workflow'
include { FUSIONREPORT_DETECT           }   from '../modules/nf-core/fusionreport/detect/main'
include { FASTQC                        }   from '../modules/nf-core/fastqc/main'
include { MULTIQC                       }   from '../modules/nf-core/multiqc/main'
include { STAR_ALIGN                    }   from '../modules/nf-core/star/align/main'
include { SALMON_QUANT                  }   from '../modules/nf-core/salmon/quant/main'
include { SAMTOOLS_CONVERT              }   from '../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_INDEX                }   from '../modules/nf-core/samtools/index/main'
include { paramsSummaryMap              }   from 'plugin/nf-schema'
include { FASTQ_ALIGN_STAR              }   from '../subworkflows/nf-core/fastq_align_star'
include { paramsSummaryMultiqc          }   from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML        }   from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText        }   from '../subworkflows/local/utils_nfcore_rnafusion_pipeline'
include { validateInputSamplesheet      }   from '../subworkflows/local/utils_nfcore_rnafusion_pipeline'
include { ARRIBA_ARRIBA                 }   from '../modules/nf-core/arriba/arriba/main'
include { CTATSPLICING_STARTOCANCERINTRONS } from '../modules/nf-core/ctatsplicing/startocancerintrons'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNAFUSION {


    take:
    ch_samplesheet_input    // channel: samplesheet read in from --input
    tools                   // list: a list of tools to run

    main:

    def ch_versions = Channel.empty()
    def ch_multiqc_files = Channel.empty()

    //
    // Create references if necessary
    //

    BUILD_REFERENCES(tools)
    ch_versions = ch_versions.mix(BUILD_REFERENCES.out.versions)

    if (!params.references_only) {

        def ch_input = ch_samplesheet_input.map { meta, fastqs, bam, bai, cram, crai, junctions, splice_junctions ->
            def align = false
            // Check if we need split junctions
            if (tools.contains("ctatsplicing") && !splice_junctions) {
                align = true
            }

            // Check if we need junctions
            if (tools.intersect(["starfusion", "ctatsplicing"]) && !junctions) {
                align = true
            }

            // Check if we need BAM or CRAM files
            if (tools.intersect(["ctatsplicing", "arriba", "stringtie", "fusioninspector"]) && !bam && !cram) {
                align = true
            }

            // Check if there are fastqs when we need to align
            if (align && !fastqs) {
                error("No fastq files found for ${meta.id}. Either provide fastq files to align or provide a BAM/CRAM file, a junctions file and a split junctions file.")
            }
            def new_meta = meta + [align:align, seq_center:meta.seq_center ?: params.seq_center, seq_platform:meta.seq_platform ?: params.seq_platform]
            return [ new_meta, fastqs, bam, bai, cram, crai, junctions, splice_junctions ]
        }
        .tap { ch_samplesheet }
        .multiMap { meta, fastqs, bam, bai, cram, crai, junctions, splice_junctions ->
            fastqs:             [ meta, fastqs ]
            bam:                [ meta, bam, bai ]
            cram:               [ meta, cram, crai ]
            junctions:          [ meta, junctions ]
            splice_junctions:   [ meta, splice_junctions ]
        }

        // Define which fastqs need to be processes (all analysis that's not aligning)
        def fastq_tools = ["salmon", "fusioninspector", "fusioncatcher"]
        selected_fastq_tools = tools.intersect(fastq_tools)
        def ch_fastqs_to_process = ch_input.fastqs.branch { meta, fastqs ->
            if (!fastqs && selected_fastq_tools) {
                log.warn("Fastq files not found for sample '${meta.id}'. Skipping the following tools for this sample: ${selected_fastq_tools.join(', ')}")
            }
            found: fastqs
            not_found: !fastqs
        }

        // Convert CRAM to BAM when needed (when tools that don't support CRAM are used and when the sample isn't aligned)
        def only_bam_tools = ["ctatsplicing", "stringtie", "fusioninspector"]
        def ch_aligned_inputs = ch_input.bam.filter { meta, file, _bai -> file && !meta.align }
        if(tools.intersect(only_bam_tools)) {
            SAMTOOLS_CONVERT(
                ch_input.cram.filter { meta, file, _crai -> file && !meta.align },
                BUILD_REFERENCES.out.fasta,
                BUILD_REFERENCES.out.fai
            )
            ch_aligned_inputs = ch_aligned_inputs.mix(
                SAMTOOLS_CONVERT.out.bam.join(SAMTOOLS_CONVERT.out.bai, failOnMismatch:true, failOnDuplicate:true)
            )
        } else {
            ch_aligned_inputs = ch_aligned_inputs.mix(ch_input.cram.filter { meta, file, _crai -> file && !meta.align })
        }

        //
        // SUBWORKFLOW: Read QC and trimming (nf-core)
        //

        def ch_reads = Channel.empty()

        // adapter_fasta as a channel (or empty if not provided)
        def ch_adapter_fasta = params.adapter_fasta ? Channel.fromPath(params.adapter_fasta).collect() : []

        // disable umi usage in this subworkflow for this pipeline
        def with_umi         = false
        def skip_umi_extract = false
        def umi_discard_read = 0

        // if 'fastp' isn't selected, we still run the subworkflow but skip trimming
        def skip_trimming    = (!tools.contains("fastp"))

        // optional fastp output controls + minimum reads after trimming
        def save_trimmed_fail = params.save_trimmed_fail ?: false
        def save_merged       = params.save_merged ?: false
        def min_trimmed_reads = (params.min_trimmed_reads ?: 1) as Integer

        FASTQ_FASTQC_UMITOOLS_FASTP(
            ch_fastqs_to_process.found,  // reads: [ val(meta), [fastqs] ]
            params.skip_qc,              // skip_fastqc
            with_umi,                    // with_umi
            skip_umi_extract,            // skip_umi_extract
            umi_discard_read,            // umi_discard_read (0,1,2)
            skip_trimming,               // skip_trimming
            ch_adapter_fasta,            // adapter_fasta
            save_trimmed_fail,           // save_trimmed_fail
            save_merged,                 // save_merged
            min_trimmed_reads            // min_trimmed_reads
        )

        ch_reads    = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads
        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.versions)

        ch_sbwf_fastp_mqc = Channel.empty()
            .mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_raw_zip.map { it[1] })
            .mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_trim_zip.map { it[1] })
            .mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_html.map { it[1] })
            .mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_json.map { it[1] })
            .ifEmpty([])

        ch_multiqc_files = ch_multiqc_files.mix(ch_sbwf_fastp_mqc)

        //
        // MODULE: SALMON_QUANT
        //

        if(tools.contains("salmon") && !params.skip_qc) {
            SALMON_QUANT(
                ch_reads,
                BUILD_REFERENCES.out.salmon_index,
                BUILD_REFERENCES.out.gtf.map{ it -> it[1] },
                [],
                false,
                'A'
            )
            ch_multiqc_files = ch_multiqc_files.mix(SALMON_QUANT.out.json_info.collect{it[1]})
            ch_versions      = ch_versions.mix(SALMON_QUANT.out.versions)
        }

        //
        // SUBWORKFLOW: Run STAR alignment
        //

        // Define which fastqs need to be aligned
        def ch_fastqs_to_align = ch_reads
            .filter { meta, _fastqs -> meta.align }

        // Add the alignment files to the correct channel if their fastqs aren't aligned
        def ch_aligned_reads        = ch_aligned_inputs
        def ch_star_junctions       = ch_input.junctions.filter { meta, file -> file && !meta.align }
        def ch_star_splice_junctions = ch_input.splice_junctions.filter { meta, file -> file && !meta.align }
        if(tools.intersect(["ctatsplicing", "arriba", "starfusion", "stringtie"])) {
            FASTQ_ALIGN_STAR(
                ch_fastqs_to_align,
                BUILD_REFERENCES.out.starindex_ref,
                BUILD_REFERENCES.out.gtf,
                params.star_ignore_sjdbgtf,
                ch_fastqs_to_align.map { meta, _fastqs -> meta.seq_platform },
                ch_fastqs_to_align.map { meta, _fastqs -> meta.seq_center },
                BUILD_REFERENCES.out.fasta,
                [[:], []]
            )
            SAMTOOLS_INDEX(FASTQ_ALIGN_STAR.out.bam_sorted_aligned)
            ch_bam_bai = FASTQ_ALIGN_STAR.out.bam_sorted_aligned
                .join(SAMTOOLS_INDEX.out.bai, failOnMismatch:true, failOnDuplicate:true)
            ch_versions             = ch_versions.mix(FASTQ_ALIGN_STAR.out.versions)
            ch_aligned_reads        = ch_aligned_reads.mix(ch_bam_bai)
            ch_star_junctions       = ch_star_junctions.mix(FASTQ_ALIGN_STAR.out.junctions)
            ch_star_splice_junctions = ch_star_splice_junctions.mix(FASTQ_ALIGN_STAR.out.spl_junc_tabs)
            ch_multiqc_files        = ch_multiqc_files.mix(FASTQ_ALIGN_STAR.out.log_final.collect{it[1]}.ifEmpty([]))
            ch_multiqc_files        = ch_multiqc_files.mix(FASTQ_ALIGN_STAR.out.gene_count.collect{it[1]}.ifEmpty([]))
        }

        //
        // MODULE: Run CTAT-SPLICING
        //

        if(tools.contains("ctatsplicing")) {
            def ch_ctatsplicing_input = ch_star_splice_junctions
                .join(ch_star_junctions, failOnMismatch:true, failOnDuplicate:true)
                .join(ch_aligned_reads, failOnMismatch:true, failOnDuplicate:true)
                .map { meta, split_junction, junction, bam, bai ->
                    [ meta, split_junction, junction, bam, bai ]
                }

            CTATSPLICING_STARTOCANCERINTRONS(
                ch_ctatsplicing_input,
                BUILD_REFERENCES.out.starfusion_ref
            )

            ch_versions = ch_versions.mix(CTATSPLICING_STARTOCANCERINTRONS.out.versions.first())
        }

        //
        // MODULE: Run Arriba
        //

        // TODO: improve how params.arriba_fusions would avoid running arriba module. Maybe imputed from samplesheet?

        def fusions_created = false
        def ch_arriba_fusions = ch_samplesheet.map { it -> [it[0], []] } // Set arriba fusions to empty by default

        if (tools.contains("arriba")) {
            fusions_created = true

            if (params.arriba_fusions) {
                ch_arriba_fusions = ch_aligned_reads.map { meta, _bam, _bai -> meta }
                    .combine(Channel.value(file(params.arriba_fusions, checkIfExists: true)))
                    .map { meta, fusion_file -> [meta, fusion_file] }
            } else {
                ARRIBA_ARRIBA(
                    ch_aligned_reads.map { meta, bam, _bai -> [meta, bam] },
                    BUILD_REFERENCES.out.fasta,
                    BUILD_REFERENCES.out.gtf,
                    BUILD_REFERENCES.out.arriba_ref_blacklist,
                    BUILD_REFERENCES.out.arriba_ref_known_fusions,
                    BUILD_REFERENCES.out.arriba_ref_cytobands,
                    BUILD_REFERENCES.out.arriba_ref_protein_domains
                )
                ch_arriba_fusions = ARRIBA_ARRIBA.out.fusions
                ch_versions = ch_versions.mix(ARRIBA_ARRIBA.out.versions)
            }
        }


        //
        // MODULE: Run StarFusion
        //

        def ch_starfusion_fusions = ch_samplesheet.map { it -> [it[0], []] } // Set starfusion fusions to empty by default
        if (tools.contains("starfusion")) {
            fusions_created = true

            if (params.starfusion_fusions) {
                def fusions = file(params.starfusion_fusions, checkIfExists:true)
                ch_starfusion_fusions = ch_star_junctions.map { meta, _junc -> [ meta, fusions ] }
            } else {
                STARFUSION_DETECT(
                    ch_star_junctions.map { meta, junc -> [ meta, [], junc ] },
                    BUILD_REFERENCES.out.starfusion_ref.map { it -> it[1] }
                )
                ch_versions = ch_versions.mix(STARFUSION_DETECT.out.versions)
                ch_starfusion_fusions = STARFUSION_DETECT.out.fusions
            }
        }


        //
        // SUBWORKFLOW: Run FusionCatcher
        //

        def ch_fusioncatcher_fusions = ch_samplesheet.map { it -> [it[0], []] } // Set fusioncatcher fusions to empty by default
        if(tools.contains("fusioncatcher")) {
            fusions_created = true
        fusioncatcher_trimming = params.trim_tail_fusioncatcher != 0
        FUSIONCATCHER_WORKFLOW (
            ch_reads,
            fusioncatcher_trimming,
            params.adapter_fasta,
            BUILD_REFERENCES.out.fusioncatcher_ref,       // channel [ meta, path       ]
            params.fusioncatcher_fusions
        )
            ch_versions = ch_versions.mix(FUSIONCATCHER_WORKFLOW.out.versions)
            // Add output of fusioncatcher to a channel + add empty entries for the samples that could not be run
            ch_fusioncatcher_fusions = FUSIONCATCHER_WORKFLOW.out.fusions.mix(ch_fastqs_to_process.not_found)
        }

        //
        // SUBWORKFLOW: Run Stringtie
        //

        if(tools.contains("stringtie")) {
            STRINGTIE_WORKFLOW (
                ch_aligned_reads.map { meta, bam, _bai -> [meta, bam]},
                BUILD_REFERENCES.out.gtf
            )
            ch_versions = ch_versions.mix(STRINGTIE_WORKFLOW.out.versions)
        }

        //
        // SUBWORKFLOW: Run FusionReport
        //

        def ch_fusion_list = Channel.empty()
        def ch_fusion_list_filtered = Channel.empty()
        def ch_fusionreport_report = Channel.empty()
        def ch_fusionreport_csv = Channel.empty()
        if (!params.skip_vis && tools.contains("fusionreport")) {
            if (!fusions_created) {
                error("Could not find any fusion files. Please generate some with `--tools arriba`, `--tools starfusion` and/or `--tools fusioncatcher`")
            }

            def ch_fusions = ch_arriba_fusions
                .join(ch_starfusion_fusions, failOnMismatch:true, failOnDuplicate:true)
                .join(ch_fusioncatcher_fusions, failOnMismatch:true, failOnDuplicate:true)

            FUSIONREPORT_DETECT(
                ch_fusions,
                BUILD_REFERENCES.out.fusionreport_ref,
                params.tools_cutoff
            )

            ch_versions             = ch_versions.mix(FUSIONREPORT_DETECT.out.versions)
            ch_fusion_list          = FUSIONREPORT_DETECT.out.fusion_list
            ch_fusion_list_filtered = FUSIONREPORT_DETECT.out.fusion_list_filtered
            ch_fusionreport_report  = FUSIONREPORT_DETECT.out.report
            ch_fusionreport_csv     = FUSIONREPORT_DETECT.out.csv
        } else if(params.fusioninspector_fusions) {
            def input_fusions       = file(params.fusioninspector_fusions, checkIfExists:true)
            ch_fusion_list          = ch_reads.map { it -> [ it[0], input_fusions ] }
            ch_fusion_list_filtered = ch_fusion_list
            ch_fusionreport_csv     = null
            ch_fusionreport_report  = null
        } else if(tools.contains("fusioninspector")) {
            error("Could not find any valid fusions for fusioninspector input. Please provide some via --fusioninspector_fusions or generate them with `--tools arriba`, `--tools starfusion` and/or `--tools fusioncatcher` with --skip_vis disabled and `--tools fusionreport enabled")
        }

        //
        // SUBWORKFLOW: Run FusionInspector
        //

        if (tools.contains("fusioninspector")) {
            FUSIONINSPECTOR_WORKFLOW (
                ch_reads,
                ch_fusion_list,
                ch_fusion_list_filtered,
                ch_fusionreport_report,
                ch_fusionreport_csv,
                ch_aligned_reads,
                BUILD_REFERENCES.out.gtf,
                BUILD_REFERENCES.out.arriba_ref_protein_domains,
                BUILD_REFERENCES.out.arriba_ref_cytobands,
                BUILD_REFERENCES.out.hgnc_ref,
                BUILD_REFERENCES.out.hgnc_date,
                BUILD_REFERENCES.out.starfusion_ref,
                params.skip_vis,
                params.skip_vcf,
                params.tools_cutoff,
                params.whitelist
            )
            ch_versions      = ch_versions.mix(FUSIONINSPECTOR_WORKFLOW.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(FUSIONINSPECTOR_WORKFLOW.out.ch_arriba_visualisation.collect{it[1]}.ifEmpty([]))
        }

        //
        // SUBWORKFLOW: Run QC
        //

        if(!params.skip_qc) {
            QC_WORKFLOW (
                ch_aligned_reads.map { meta, bam, _bai -> [meta, bam] },
                BUILD_REFERENCES.out.refflat,
                BUILD_REFERENCES.out.fasta,
                BUILD_REFERENCES.out.fai,
                BUILD_REFERENCES.out.rrna_interval
            )
            ch_versions      = ch_versions.mix(QC_WORKFLOW.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(QC_WORKFLOW.out.rnaseq_metrics.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(QC_WORKFLOW.out.duplicate_metrics.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(QC_WORKFLOW.out.insertsize_metrics.collect{it[1]})
        }
    }
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'rnafusion_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //

    def ch_multiqc_output = Channel.empty()
    if(!params.skip_qc && !params.references_only) {
        ch_multiqc_config        = Channel.fromPath(
            "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config = params.multiqc_config ?
            Channel.fromPath(params.multiqc_config, checkIfExists: true) :
            Channel.empty()
        ch_multiqc_logo          = params.multiqc_logo ?
            Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
            Channel.empty()

        summary_params      = paramsSummaryMap(
            workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
        ch_multiqc_files = ch_multiqc_files.mix(
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
            file(params.multiqc_methods_description, checkIfExists: true) :
            file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        ch_methods_description                = Channel.value(
            methodsDescriptionText(ch_multiqc_custom_methods_description))

        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
        ch_multiqc_files = ch_multiqc_files.mix(
            ch_methods_description.collectFile(
                name: 'methods_description_mqc.yaml',
                sort: true
            )
        )

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            [],
            []
        )
        ch_multiqc_output = MULTIQC.out.report.toList()
    }


    emit:
    multiqc_report =  ch_multiqc_output // channel: /path/to/multiqc_report.html
    versions       = ch_versions        // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
