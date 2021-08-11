/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowRnafusion.initialise(params, log)

def checkPathParamList = [
    params.input,
    params.multiqc_config,
    params.fasta
]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Check alignment parameters
def prepareToolIndices  = []
prepareToolIndices << params.aligner

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()
def publish_genome_options = params.save_reference ? [publish_dir: 'genome']       : [publish_files: false]
def publish_index_options  = params.save_reference ? [publish_dir: 'genome/index'] : [publish_files: false]
def star_genomegenerate_options = modules['star_genomegenerate']
def starfusion_genome_options = modules['starfusion_download']
def fusioncatcher_genome_options = modules['fusioncatcher_download']
if (!params.save_reference)     { star_genomegenerate_options['publish_files'] = false }

// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )
//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'                addParams( options: [:] )
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'          addParams( genome_options: publish_genome_options, index_options: publish_index_options, star_index_options: star_genomegenerate_options, starfusion_untar_options: starfusion_genome_options, starfusion_download_options: starfusion_genome_options, fusioncatcher_download_options: fusioncatcher_genome_options)
/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                }   from '../modules/nf-core/modules/fastqc/main'           addParams( options: modules['fastqc'] )
include { MULTIQC               }   from '../modules/nf-core/modules/multiqc/main'          addParams( options: multiqc_options   )
include { FUSION_STAR_ARRIBA    }   from '../subworkflows/nf-core/fusion_star_arriba'       addParams( star_align_options: modules['star_align'], arriba_options: modules['arriba_fusion'])
include { STARFUSION            }   from '../modules/local/starfusion/detection/main'       addParams( options: modules['starfusion'] )
include { FUSIONCATCHER         }   from '../modules/local/fusioncatcher/detection/main'    addParams( options: modules['fusioncatcher'] )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow RNAFUSION {

    // SUBWORKFLOW: Uncompress and prepare reference genome files
    PREPARE_GENOME (
        prepareToolIndices
    )

    ch_software_versions = Channel.empty()
    ch_reads = Channel.empty()
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_reads = INPUT_CHECK.out.reads
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_reads
    )
    ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))

    //
    // MODULE: Pipeline reporting
    //
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )

    // MODULE: MultiQC

    workflow_summary    = WorkflowRnafusion.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))

    // Run STAR alignment and Arriba
    if (params.arriba){
        ch_genome_bam     = Channel.empty()
        ch_arriba_fusions = Channel.empty()
        FUSION_STAR_ARRIBA (
            ch_reads,
            PREPARE_GENOME.out.star_index,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gtf
        )
        ch_genome_bam     = FUSION_STAR_ARRIBA.out.bam
        ch_arriba_fusions = FUSION_STAR_ARRIBA.out.fusions
    }

    //Run STAR fusion
    if (params.starfusion){
        ch_star_fusions = Channel.empty()
        STARFUSION (
            ch_reads,
            PREPARE_GENOME.out.starfusion_resource
        )
        ch_starfusion_fusions = STARFUSION.out.fusions
    }

    //Run FusionCatcher
    if (params.fusioncatcher){
        ch_fusioncather_fusions = Channel.empty()
        FUSIONCATCHER (
            ch_reads,
            PREPARE_GENOME.out.fusioncatcher_resource
        )
        ch_fusioncatcher_fusions = FUSIONCATCHER.out.fusions
    }
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
