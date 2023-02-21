/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowRnafusion.initialise(params, log)

// Check mandatory parameters

if (file(params.input).exists() || params.build_references) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet does not exist or was not specified!' }
if (params.fusioninspector_only && !params.fusioninspector_fusions) { exit 1, 'Parameter --fusioninspector_fusions PATH_TO_FUSION_LIST expected with parameter --fusioninspector_only'}

ch_chrgtf = params.starfusion_build ? file(params.chrgtf) : file("${params.starfusion_ref}/ref_annot.gtf")
ch_starindex_ref = params.starfusion_build ? params.starindex_ref : "${params.starfusion_ref}/ref_genome.fa.star.idx"
ch_starindex_ensembl_ref = params.starindex_ref
ch_refflat = params.starfusion_build ? file(params.refflat) : "${params.ensembl_ref}/ref_annot.gtf.refflat"
ch_rrna_interval = params.starfusion_build ? file(params.rrna_intervals) : "${params.ensembl_ref}/ref_annot.interval_list"


def checkPathParamList = [
    params.fasta,
    params.fai,
    params.gtf,
    ch_chrgtf,
    params.transcript,
    ch_refflat,
    ch_rrna_interval
]


for (param in checkPathParamList) if ((param) && !params.build_references) file(param, checkIfExists: true)
if (params.fasta[0,1] == "s3") {
    log.info "INFO: s3 path detected, check for absolute path and trailing '/' not performed"
}
else {
    for (param in checkPathParamList) if ((param.toString())!= file(param).toString() && !params.build_references) { exit 1, "Problem with ${param}: ABSOLUTE PATHS are required! Check for trailing '/' at the end of paths too." }
}
if ((params.squid || params.all) && params.ensembl_version == 105) { exit 1, 'Ensembl version 105 is not supported by squid' }

ch_fasta = file(params.fasta)
ch_gtf = file(params.gtf)
ch_transcript = file(params.transcript)
ch_fai = file(params.fai)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { INPUT_CHECK                   }   from '../subworkflows/local/input_check'
include { TRIM_WORKFLOW                 }   from '../subworkflows/local/trim_workflow'
include { ARRIBA_WORKFLOW               }   from '../subworkflows/local/arriba_workflow'
include { PIZZLY_WORKFLOW               }   from '../subworkflows/local/pizzly_workflow'
include { QC_WORKFLOW                   }   from '../subworkflows/local/qc_workflow'
include { SQUID_WORKFLOW                }   from '../subworkflows/local/squid_workflow'
include { STARFUSION_WORKFLOW           }   from '../subworkflows/local/starfusion_workflow'
include { STRINGTIE_WORKFLOW            }   from '../subworkflows/local/stringtie_workflow'
include { FUSIONCATCHER_WORKFLOW        }   from '../subworkflows/local/fusioncatcher_workflow'
include { FUSIONINSPECTOR_WORKFLOW      }   from '../subworkflows/local/fusioninspector_workflow'
include { FUSIONREPORT_WORKFLOW         }   from '../subworkflows/local/fusionreport_workflow'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq/main'
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow RNAFUSION {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    .reads
    .map {
        meta, fastq ->
            def meta_clone = meta.clone()
            meta_clone.id = meta_clone.id.split('_')[0..-2].join('_')
            [ meta_clone, fastq ] }
    .groupTuple(by: [0])
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))


    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_cat_fastq
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    TRIM_WORKFLOW (
        ch_cat_fastq
    )
    ch_reads_fusioncatcher = TRIM_WORKFLOW.out.ch_reads_fusioncatcher
    ch_reads_all = TRIM_WORKFLOW.out.ch_reads_all


    // Run STAR alignment and Arriba
    ARRIBA_WORKFLOW (
        ch_reads_all,
        ch_gtf,
        ch_fasta,
        ch_starindex_ensembl_ref
    )
    ch_versions = ch_versions.mix(ARRIBA_WORKFLOW.out.versions.first().ifEmpty(null))

    // Run pizzly/kallisto

    PIZZLY_WORKFLOW (
        ch_reads_all,
        ch_gtf,
        ch_transcript
    )
    ch_versions = ch_versions.mix(PIZZLY_WORKFLOW.out.versions.first().ifEmpty(null))


// Run squid

    SQUID_WORKFLOW (
        ch_reads_all,
        ch_gtf,
        ch_starindex_ensembl_ref,
        ch_fasta
    )
    ch_versions = ch_versions.mix(SQUID_WORKFLOW.out.versions.first().ifEmpty(null))


//Run STAR fusion
    STARFUSION_WORKFLOW (
        ch_reads_all,
        ch_chrgtf,
        ch_starindex_ref
    )
    ch_versions = ch_versions.mix(STARFUSION_WORKFLOW.out.versions.first().ifEmpty(null))


//Run fusioncatcher
    FUSIONCATCHER_WORKFLOW (
        ch_reads_fusioncatcher
    )
    ch_versions = ch_versions.mix(FUSIONCATCHER_WORKFLOW.out.versions.first().ifEmpty(null))


//Run stringtie
    STRINGTIE_WORKFLOW (
        STARFUSION_WORKFLOW.out.bam_sorted,
        ch_chrgtf
    )
    ch_versions = ch_versions.mix(STRINGTIE_WORKFLOW.out.versions.first().ifEmpty(null))


    //Run fusion-report
    FUSIONREPORT_WORKFLOW (
        ch_reads_all,
        params.fusionreport_ref,
        ARRIBA_WORKFLOW.out.fusions,
        PIZZLY_WORKFLOW.out.fusions,
        SQUID_WORKFLOW.out.fusions,
        STARFUSION_WORKFLOW.out.fusions,
        FUSIONCATCHER_WORKFLOW.out.fusions
    )
    ch_versions = ch_versions.mix(FUSIONREPORT_WORKFLOW.out.versions.first().ifEmpty(null))


    //Run fusionInpector
    FUSIONINSPECTOR_WORKFLOW (
        ch_reads_all,
        FUSIONREPORT_WORKFLOW.out.fusion_list,
        FUSIONREPORT_WORKFLOW.out.fusion_list_filtered
    )
    ch_versions = ch_versions.mix(FUSIONINSPECTOR_WORKFLOW.out.versions.first().ifEmpty(null))


    //QC
    QC_WORKFLOW (
        STARFUSION_WORKFLOW.out.bam_sorted,
        ch_chrgtf,
        ch_refflat,
        ch_fasta,
        ch_fai,
        ch_rrna_interval
    )
    ch_versions = ch_versions.mix(QC_WORKFLOW.out.versions.first().ifEmpty(null))

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )


    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowRnafusion.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowRnafusion.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_WORKFLOW.out.qualimap_qc.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(STARFUSION_WORKFLOW.out.star_stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_WORKFLOW.out.rnaseq_metrics.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_WORKFLOW.out.duplicate_metrics.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    multiqc_report       = MULTIQC.out.report.toList()
    ch_versions          = ch_versions.mix(MULTIQC.out.versions)
}




/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
