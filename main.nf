#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/rnafusion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/rnafusion
    Website: https://nf-co.re/rnafusion
    Slack  : https://nfcore.slack.com/channels/rnafusion
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_rnafusion_pipeline'

params.fasta                        = getGenomeAttribute("fasta")
params.fai                          = getGenomeAttribute("fai")
params.gtf                          = getGenomeAttribute("gtf")
params.refflat                      = getGenomeAttribute("refflat")
params.rrna_intervals               = getGenomeAttribute("rrna_intervals")
params.arriba_ref_blacklist         = getGenomeAttribute("arriba_ref_blacklist")
params.arriba_ref_cytobands         = getGenomeAttribute("arriba_ref_cytobands")
params.arriba_ref_known_fusions     = getGenomeAttribute("arriba_ref_known_fusions")
params.arriba_ref_protein_domains   = getGenomeAttribute("arriba_ref_protein_domains")
params.fusioncatcher_ref            = getGenomeAttribute("fusioncatcher_ref")
params.hgnc_ref                     = getGenomeAttribute("hgnc_ref")
params.hgnc_date                    = getGenomeAttribute("hgnc_date")
params.salmon_index                 = getGenomeAttribute("salmon_index")
params.starfusion_ref               = getGenomeAttribute("starfusion_ref")
params.starindex_ref                = getGenomeAttribute("starindex_ref")
params.fusionreport_ref             = getGenomeAttribute("fusionreport_ref")

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_rnafusion_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_rnafusion_pipeline'
include { RNAFUSION               } from './workflows/rnafusion'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        args,
        params.outdir,
    )

    def tools = params.tools.tokenize(",")
    if (tools.contains("all")) {
        def json = new groovy.json.JsonSlurper().parseText(file("${projectDir}/nextflow_schema.json").text)
        def pattern = json.get('$defs')?.get('input_output_options')?.get('properties')?.get('tools')?.get('pattern')
        if (!pattern) {
            error("Could not fetch the allowed tools from the JSON schema, please check the code. If you see this as a pipeline user, please contact the developers instead.")
        }
        tools = pattern.replace('^((', "").replace(')?,?)*(?<!,)$', "").tokenize("|") - "all"
    }
    log.debug("Rnafusion tools to run: ${tools}")

    def profiles = workflow.profile
    if ((profiles.contains("conda") || profiles.contains("mamba")) && (tools.contains("ctatsplicing"))) {
        error("Conda or Mamba runs are not supported when ctatsplicing is in `--tools`")
    }

    if (tools.contains("fusioncatcher") && (!params.fusioncatcher_ref || !file(params.fusioncatcher_ref).exists())) {
        error("You have selected `fusioncatcher` in `--tools`, but did not provide an existing path to the fusioncatcher reference files with `--fusioncatcher_ref`.")
    }

    if (tools.contains("arriba")) {
        if (!params.arriba_ref_blacklist || !file(params.arriba_ref_blacklist).exists()) {
            error("You have selected `arriba` in `--tools`, but did not provide an existing path to the arriba reference blacklist file with `--arriba_ref_blacklist`.")
        }
        if (!params.arriba_ref_cytobands || !file(params.arriba_ref_cytobands).exists()) {
            error("You have selected `arriba` in `--tools`, but did not provide an existing path to the arriba reference cytobands file with `--arriba_ref_cytobands`.")
        }
        if (!params.arriba_ref_known_fusions || !file(params.arriba_ref_known_fusions).exists()) {
            error("You have selected `arriba` in `--tools`, but did not provide an existing path to the arriba reference known fusions file with `--arriba_ref_known_fusions`.")
        }
        if (!params.arriba_ref_protein_domains || !file(params.arriba_ref_protein_domains).exists()) {
            error("You have selected `arriba` in `--tools`, but did not provide an existing path to the arriba reference protein domains file with `--arriba_ref_protein_domains`.")
        }
    }

    //
    // WORKFLOW: Run main workflow
    //
    RNAFUSION(
        PIPELINE_INITIALISATION.out.samplesheet,
        tools
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        RNAFUSION.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
