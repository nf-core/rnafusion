//
// Subworkflow with functionality specific to the nf-core/rnafusion pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //

    Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map { meta, fastq_1, fastq_2, bam, bai, cram, crai, junctions, splice_junctions, strandedness ->
            def meta_fastqs = []
            if (!fastq_1) {
                meta_fastqs = [ meta, [] ]
            } else if (!fastq_2) {
                meta_fastqs = [ meta + [single_end:true], [ fastq_1 ] ]
            } else {
                meta_fastqs = [ meta + [single_end:false], [ fastq_1, fastq_2 ] ]
            }
            return [ meta.id ] + meta_fastqs + [ bam, bai, cram, crai, junctions, splice_junctions, strandedness ]
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .map {
            meta, fastqs, bam, bai, cram, crai, junctions, splice_junctions ->
                return [ meta, fastqs.flatten(), bam, bai, cram, crai, junctions, splice_junctions ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal()
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()

    if (params.no_cosmic) {
        log.warn("Skipping COSMIC DB download from `FUSIONREPORT_DOWNLOAD` and skip using it in `FUSIONREPORT`")
    }

    def dfamParams = ['dfam_hmm', 'dfam_h3p', 'dfam_h3m', 'dfam_h3i', 'dfam_h3f']

    if (params.dfam_version) {
        def dfamPattern = "https://www.dfam.org/releases/Dfam_${params.dfam_version}/infrastructure/dfamscan/${params.species}_dfam"

        def setDfamParams = dfamParams.findAll { params[it] }

        if (setDfamParams) {
            def customParams = setDfamParams.findAll { paramName ->
                !params[paramName]?.startsWith(dfamPattern)
            }
            if (customParams) {
                def paramDetails = customParams.collect { paramName ->
                    "   --${paramName}: ${params[paramName]}"
                }.join('\n')
                def dfam_warn = "Both custom DFAM paths as well as `--dfam_version` (${params.dfam_version}) and `--species` (${params.species}) were provided.\n" +
                    "Custom DFAM paths parameters provided:\n${paramDetails}\n" +
                    "The pipeline will prioritize these custom files specified with `--${customParams}` and **will NOT** construct these URLs based on `--dfam_version` nor `--species`.\n" +
                    "   - If you intend to use custom DFAM files, please ensure that all `--dfam_h*` parameters point to full and valid paths.\n" +
                    "   - If you prefer to let the pipeline build the DFAM URLs automatically, omit `--dfam_h*` and instead provide only `--dfam_version` and `--species`."
                log.warn(dfam_warn)
            }
        }
    }

    if (params.pfam_version){
        def pfamPattern = "http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam${params.pfam_version}/Pfam-A"

        if (!(params.pfam_file?.startsWith(pfamPattern))) {
            def pfam_warn = "Both `--pfam_file` (${params.pfam_file}) and `--pfam_version` (${params.pfam_version}) were provided.\n" +
                    "The pipeline will prioritize the custom file from `--pfam_file` and **will NOT** construct the URL based on `--pfam_version`.\n" +
                    "   - If you intend to use a custom PFAM file, please ensure that `--pfam_file` points to a full and valid path.\n" +
                    "   - If you prefer to let the pipeline build the PFAM URL automatically, omit `--pfam_file` and instead provide only `--pfam_version`."
            log.warn(pfam_warn)
        }
    }
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs, bam, bai, cram, crai, junctions, splice_junctions) = input[1..8]

    def bam_list = bam.findAll { it -> it != [] }
    def cram_list = cram.findAll { it -> it != [] }
    def junctions_list = junctions.findAll { it -> it != [] }
    def splice_junctions_list = splice_junctions.findAll { it -> it != [] }
    // Check alignment and junction files (input is a list)
    if (bam_list.size() > 1 || cram_list.size() > 1 || junctions_list.size() > 1 || splice_junctions_list.size() > 1) {
        error("Please check input samplesheet -> Only one BAM or CRAM, junctions and split junctions file is allowed per sample: ${metas[0].id}")
    }

    bam = bam_list.size() > 0 ? bam_list[0] : []
    cram = cram_list.size() > 0 ? cram_list[0] : []
    junctions = junctions_list.size() > 0 ? junctions_list[0] : []
    splice_junctions = splice_junctions_list.size() > 0 ? splice_junctions_list[0] : []

    if (bam != [] && cram != []) {
        error("Please check input samplesheet -> Using both BAM and CRAM files isn't allowed: ${metas[0].id}")
    }

    // Check that multiple runs of the same sample are of the same strandedness
    def strandedness_ok = metas.collect{ it.strandedness }.unique().size == 1
    if (!strandedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must have the same strandedness!: ${metas[0].id}")
    }

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok && fastqs) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs, bam, bai.find { it -> it != [] } ?: [], cram, crai.find { it -> it != [] } ?: [], junctions, splice_junctions ]
}

//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()

    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
