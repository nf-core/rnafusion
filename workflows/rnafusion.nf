/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { SALMON_QUANT           } from '../modules/nf-core/salmon/quant/main'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_rnafusion_pipeline'
include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { BUILD_REFERENCES              }   from '../subworkflows/local/build_references'
include { TRIM_WORKFLOW                 }   from '../subworkflows/local/trim_workflow'
include { ARRIBA_WORKFLOW               }   from '../subworkflows/local/arriba_workflow'
include { QC_WORKFLOW                   }   from '../subworkflows/local/qc_workflow'
include { STARFUSION_WORKFLOW           }   from '../subworkflows/local/starfusion_workflow'
include { STRINGTIE_WORKFLOW            }   from '../subworkflows/local/stringtie_workflow'
include { FUSIONCATCHER_WORKFLOW        }   from '../subworkflows/local/fusioncatcher_workflow'
include { FUSIONINSPECTOR_WORKFLOW      }   from '../subworkflows/local/fusioninspector_workflow'
include { FUSIONREPORT_WORKFLOW         }   from '../subworkflows/local/fusionreport_workflow'
include { validateInputSamplesheet      }   from '../subworkflows/local/utils_nfcore_rnafusion_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNAFUSION {

    take:
    ch_samplesheet

    main:


    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()


    //
    // Create references if necessary
    //

    BUILD_REFERENCES()


    //
    // Create channel from input file provided through params.input
    //
    Channel
        .fromSamplesheet("input")
        .map {
            meta, fastq_1, fastq_2, strandedness ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map {
            validateInputSamplesheet(it)
        }
        .branch {
            meta, fastqs ->
                single  : fastqs.size() == 1
                    return [ meta, fastqs.flatten() ]
                multiple: fastqs.size() > 1
                    return [ meta, fastqs.flatten() ]
        }
        .set { ch_fastq }

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)


    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_cat_fastq
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions)


    TRIM_WORKFLOW (
        ch_cat_fastq
    )
    ch_reads_fusioncatcher = TRIM_WORKFLOW.out.ch_reads_fusioncatcher
    ch_reads_all = TRIM_WORKFLOW.out.ch_reads_all
    ch_versions = ch_versions.mix(TRIM_WORKFLOW.out.versions)


//     SALMON_QUANT( ch_reads_all, ch_salmon_index.map{ meta, index ->  index  }, ch_gtf.map{ meta, gtf ->  gtf  }, [], false, 'A')


//     //
//     // SUBWORKFLOW:  Run STAR alignment and Arriba
//     //
//     ARRIBA_WORKFLOW (
//         ch_reads_all,
//         ch_gtf,
//         ch_fasta,
//         ch_starindex_ensembl_ref,
//         ch_arriba_ref_blacklist,
//         ch_arriba_ref_known_fusions,
//         ch_arriba_ref_protein_domains
//     )
//     ch_versions = ch_versions.mix(ARRIBA_WORKFLOW.out.versions)


// //Run STAR fusion
//     STARFUSION_WORKFLOW (
//         ch_reads_all,
//         ch_chrgtf,
//         ch_starindex_ref,
//         ch_fasta
//     )
//     ch_versions = ch_versions.mix(STARFUSION_WORKFLOW.out.versions)


// //Run fusioncatcher
//     FUSIONCATCHER_WORKFLOW (
//         ch_reads_fusioncatcher
//     )
//     ch_versions = ch_versions.mix(FUSIONCATCHER_WORKFLOW.out.versions)


// //Run stringtie
//     STRINGTIE_WORKFLOW (
//         STARFUSION_WORKFLOW.out.ch_bam_sorted,
//         ch_chrgtf
//     )
//     ch_versions = ch_versions.mix(STRINGTIE_WORKFLOW.out.versions)


//     //Run fusion-report
//     FUSIONREPORT_WORKFLOW (
//         ch_reads_all,
//         ch_fusionreport_ref,
//         ARRIBA_WORKFLOW.out.fusions,
//         STARFUSION_WORKFLOW.out.fusions,
//         FUSIONCATCHER_WORKFLOW.out.fusions
//     )
//     ch_versions = ch_versions.mix(FUSIONREPORT_WORKFLOW.out.versions)


//     //Run fusionInpector
//     FUSIONINSPECTOR_WORKFLOW (
//         ch_reads_all,
//         FUSIONREPORT_WORKFLOW.out.fusion_list,
//         FUSIONREPORT_WORKFLOW.out.fusion_list_filtered,
//         FUSIONREPORT_WORKFLOW.out.report,
//         FUSIONREPORT_WORKFLOW.out.csv,
//         STARFUSION_WORKFLOW.out.ch_bam_sorted_indexed,
//         ch_chrgtf,
//         ch_arriba_ref_protein_domains,
//         ch_arriba_ref_cytobands,
//         ch_hgnc_ref,
//         ch_hgnc_date
//     )
//     ch_versions = ch_versions.mix(FUSIONINSPECTOR_WORKFLOW.out.versions)


//     //QC
//     QC_WORKFLOW (
//         ch_reads_all,
//         STARFUSION_WORKFLOW.out.ch_bam_sorted,
//         STARFUSION_WORKFLOW.out.ch_bam_sorted_indexed,
//         ch_chrgtf,
//         ch_refflat,
//         ch_fasta,
//         ch_fai,
//         ch_rrna_interval
//     )
//     ch_versions = ch_versions.mix(QC_WORKFLOW.out.versions)

//     //
//     // Collate and save software versions
//     //
//     softwareVersionsToYAML(ch_versions)
//         .collectFile(
//             storeDir: "${params.outdir}/pipeline_info",
//             name: 'nf_core_pipeline_software_mqc_versions.yml',
//             sort: true,
//             newLine: true
//         ).set { ch_collated_versions }


//     //
//     // MODULE: MultiQC
//     //
//     ch_multiqc_config        = Channel.fromPath(
//         "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
//     ch_multiqc_custom_config = params.multiqc_config ?
//         Channel.fromPath(params.multiqc_config, checkIfExists: true) :
//         Channel.empty()
//     ch_multiqc_logo          = params.multiqc_logo ?
//         Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
//         Channel.empty()

//     summary_params      = paramsSummaryMap(
//         workflow, parameters_schema: "nextflow_schema.json")
//     ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

//     ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
//         file(params.multiqc_methods_description, checkIfExists: true) :
//         file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
//     ch_methods_description                = Channel.value(
//         methodsDescriptionText(ch_multiqc_custom_methods_description))

//     ch_multiqc_files = ch_multiqc_files.mix(
//         ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
//     ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
//     ch_multiqc_files = ch_multiqc_files.mix(
//         ch_methods_description.collectFile(
//             name: 'methods_description_mqc.yaml',
//             sort: true
//         )
//     )
//     ch_multiqc_files                      = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
//     ch_multiqc_files                      = ch_multiqc_files.mix(TRIM_WORKFLOW.out.ch_fastp_html.collect{it[1]}.ifEmpty([]))
//     ch_multiqc_files                      = ch_multiqc_files.mix(TRIM_WORKFLOW.out.ch_fastp_json.collect{it[1]}.ifEmpty([]))
//     ch_multiqc_files                      = ch_multiqc_files.mix(TRIM_WORKFLOW.out.ch_fastqc_trimmed.collect{it[1]}.ifEmpty([]))
//     ch_multiqc_files                      = ch_multiqc_files.mix(STARFUSION_WORKFLOW.out.star_stats.collect{it[1]}.ifEmpty([]))
//     ch_multiqc_files                      = ch_multiqc_files.mix(STARFUSION_WORKFLOW.out.star_gene_count.collect{it[1]}.ifEmpty([]))
//     ch_multiqc_files                      = ch_multiqc_files.mix(QC_WORKFLOW.out.rnaseq_metrics.collect{it[1]}.ifEmpty([]))
//     ch_multiqc_files                      = ch_multiqc_files.mix(QC_WORKFLOW.out.duplicate_metrics.collect{it[1]}.ifEmpty([]))
//     ch_multiqc_files                      = ch_multiqc_files.mix(QC_WORKFLOW.out.insertsize_metrics.collect{it[1]}.ifEmpty([]))
//     ch_multiqc_files                      = ch_multiqc_files.mix(FUSIONINSPECTOR_WORKFLOW.out.ch_arriba_visualisation.collect{it[1]}.ifEmpty([]))

//     MULTIQC (
//         ch_multiqc_files.collect(),
//         ch_multiqc_config.toList(),
//         ch_multiqc_custom_config.toList(),
//         ch_multiqc_logo.toList()
//     )

//     emit:
//     multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
//     versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
