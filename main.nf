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

params {

    // Skip QC steps
    skip_qc: Boolean

    // Skip vcf creation step
    skip_vcf: Boolean

    // Skip visualisation steps
    skip_vis: Boolean

    // Path to samplesheet file containing information about the samples in the experiment.
    input: Path

    // The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.
    outdir: String

    // Email address for completion summary.
    email: String?

    // MultiQC report title. Printed as page header, used for filename if not otherwise specified.
    multiqc_title: String?

    // COSMIC username
    cosmic_username: String?

    // COSMIC password
    cosmic_passwd: String?

    // Path to reference folder
    genomes_base: String = "s3://nf-core-awsmegatests/rnafusion/references"

    // Don't automatically assign reference parameters to the correct references in --genomes_base
    genomes_ignore: Boolean

    // gencode version
    genome_gencode_version: String // specified in config

    // Comma-delimited list of tools to run
    tools: String

    // The length of the reads provided to the pipeline. This is used for the '--sjdbOverhang' option of STAR as read_length - 1. Providing 1 to this option will disable overhang handling.
    read_length: Integer = 100

    // Path to arriba reference blacklist
    arriba_ref_blacklist: String = getGenomeAttribute("arriba_ref_blacklist", params.genomes, params.genome)

    // Path to arriba reference cytobands
    arriba_ref_cytobands: String = getGenomeAttribute("arriba_ref_cytobands", params.genomes, params.genome)

    // Path to arriba reference known fusions
    arriba_ref_known_fusions: String = getGenomeAttribute("arriba_ref_known_fusions", params.genomes, params.genome)

    // Path to arriba reference protein domain
    arriba_ref_protein_domains: String = getGenomeAttribute("arriba_ref_protein_domains", params.genomes, params.genome)

    // Path to arriba output
    arriba_fusions: Path?

    // Path to fusioncatcher output
    fusioncatcher_fusions: Path?

    // Use limitSjdbInsertNsj with int for fusioncatcher
    fusioncatcher_limitSjdbInsertNsj: Integer = 2000000

    // Path to fusioncatcher references
    fusioncatcher_ref: String = getGenomeAttribute("fusioncatcher_ref", params.genomes, params.genome)

    // Use limitSjdbInsertNsj with int for fusioninspector STAR process
    fusioninspector_limitSjdbInsertNsj: Integer = 1000000

    // Path to a fusion list file built with format GENE1--GENE2
    fusioninspector_fusions: Path?

    // Path to fusionreport references
    fusionreport_ref: String = getGenomeAttribute("fusionreport_ref", params.genomes, params.genome)

    // Path to HGNC database file
    hgnc_ref: String = getGenomeAttribute("hgnc_ref", params.genomes, params.genome)

    // Path to HGNC timestamp file for database retrieval
    hgnc_date: String = getGenomeAttribute("hgnc_date", params.genomes, params.genome)

    // Use QIAGEN instead of SANGER to download COSMIC database
    qiagen: Boolean

    // Path to salmon index
    salmon_index: String = getGenomeAttribute("salmon_index", params.genomes, params.genome)

    // Path to starfusion output
    starfusion_fusions: Path?

    // Path to starfusion references
    starfusion_ref: String = getGenomeAttribute("starfusion_ref", params.genomes, params.genome)

    // Path to the cancer introns CSV file to create the CTAT-SPLICING reference with
    ctatsplicing_cancer_introns: Path = 'https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/CANCER_SPLICING_LIB_SUPPLEMENT/cancer_introns.GRCh38.Jun232020.tsv.gz'

    // Path to Dfam HMM database file
    dfam_hmm: Path = "https://www.dfam.org/releases/Dfam_${params.dfam_version}/infrastructure/dfamscan/${params.species}_dfam.hmm"

    // Path to Dfam H3F database file
    dfam_h3f: Path = "https://www.dfam.org/releases/Dfam_${params.dfam_version}/infrastructure/dfamscan/${params.species}_dfam.hmm.h3f"

    // Path to Dfam H3I database file
    dfam_h3i: Path = "https://www.dfam.org/releases/Dfam_${params.dfam_version}/infrastructure/dfamscan/${params.species}_dfam.hmm.h3i"

    // Path to Dfam H3M database file
    dfam_h3m: Path = "https://www.dfam.org/releases/Dfam_${params.dfam_version}/infrastructure/dfamscan/${params.species}_dfam.hmm.h3m"

    // Path to Dfam H3P database file
    dfam_h3p: Path = "https://www.dfam.org/releases/Dfam_${params.dfam_version}/infrastructure/dfamscan/${params.species}_dfam.hmm.h3p"

    // Path to Pfam database file
    pfam_file: Path = "http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam${params.pfam_version}/Pfam-A.hmm.gz"

    // URL to annotation filter rule file
    annot_filter_url: Path = 'https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/AnnotFilterRule.pm'

    // Path to starindex references
    starindex_ref: String = getGenomeAttribute("starindex_ref", params.genomes, params.genome)

    // Discard fusions identified by less than INT tools
    tools_cutoff: Integer = 1

    // Path to fusions to add to the input of fusioninspector
    whitelist: Path?

    // The amount of bases to trim at the tail of each read, none will be trimmed by default
    trim_tail: Integer = 0

    // The amount of bases to trim at the tail of each read for fusioncatcher, none will be trimmed by default
    trim_tail_fusioncatcher: Integer = 0

    // Path to adapter fasta file
    adapter_fasta: Path?

    // FASTP: Specify true to save files that failed to pass trimming thresholds
    save_trimmed_fail: Boolean

    // FASTP: Inputs with fewer than this reads will be filtered out of the "reads" output channel
    min_trimmed_reads: Integer = 1

    // FASTP: Specify true to save merged reads
    save_merged: Boolean

    // Output CRAM files instead of BAM files.
    cram: Boolean

    // Skip running the analysis, only builds the references
    references_only: Boolean

    // Path to FASTA genome file.
    fasta: String = getGenomeAttribute("fasta", params.genomes, params.genome)

    // Path to FASTA genome index file.
    fai: String = getGenomeAttribute("fai", params.genomes, params.genome)

    // Name of iGenomes reference.
    genome: String

    // Path to GTF genome file.
    gtf: String = getGenomeAttribute("gtf", params.genomes, params.genome)

    // Path to GTF genome file.
    refflat: String = getGenomeAttribute("refflat", params.genomes, params.genome)

    // Path to ribosomal interval list.
    rrna_intervals: String = getGenomeAttribute("rrna_intervals", params.genomes, params.genome)

    // Avoid using Cosmic DB (for example in clinical case applications where a paid license applies.
    no_cosmic: Boolean

    // Path to Fusion Annotation Library to be used in STARFUSION_BUILD.
    fusion_annot_lib: Path = "https://github.com/FusionAnnotator/CTAT_HumanFusionLib/releases/download/v0.3.0/fusion_lib.Mar2021.dat.gz" // path to  dat.gz CTAT genome lib // TODO: Update to latest with s3 link when available

    // Which species dfam should automatically download, default: homo_sapiens.
    species: String = 'homo_sapiens'

    // Version of dfam to use
    dfam_version: String // specified in config

    // Version of pfam to use
    pfam_version: String // specified in config

    // Git commit id for Institutional configs.
    custom_config_version: String = 'master'

    // Base directory for Institutional configs.
    custom_config_base: String = 'https://raw.githubusercontent.com/nf-core/configs/master'

    // Institutional config name.
    config_profile_name: String?

    // Institutional config description.
    config_profile_description: String?

    // Institutional config contact information.
    config_profile_contact: String?

    // Institutional config URL link.
    config_profile_url: String?

    // Display version and exit.
    version: Boolean

    // Method used to save pipeline results to output directory.
    publish_dir_mode: String = 'copy'

    // Email address for completion summary, only when pipeline fails.
    email_on_fail: String?

    // Send plain-text email instead of HTML.
    plaintext_email: Boolean

    // File size limit when attaching MultiQC reports to summary emails.
    max_multiqc_email_size: MemoryUnit = 25.MB

    // Do not use coloured log outputs.
    monochrome_logs: Boolean

    // Custom config file to supply to MultiQC.
    multiqc_config: Path?

    // Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file
    multiqc_logo: Path?

    // Custom MultiQC yaml file containing HTML including a methods description.
    multiqc_methods_description: Path?

    // Boolean whether to validate parameters against the schema at runtime
    validate_params: Boolean = true

    // Base URL or local path to location of pipeline test dataset files
    pipelines_testdata_base_path: String = 'https://raw.githubusercontent.com/nf-core/test-datasets/'

    // Suffix to add to the trace report filename. Default is the date and time in the format yyyy-MM-dd_HH-mm-ss.
    trace_report_suffix: String // specified in config

    // Display a short help message. Give a parameter name to get detailed help for that parameter.
    help = false

    // Display the full detailed help message.
    help_full: Boolean

    // Display hidden parameters in the help message (only works when --help or --help_full are provided).
    show_hidden: Boolean

    // Sequencing center, used to fill in read group CN tag in the BAM header
    seq_center: String = ""

    // Sequencing platform, used to fill in read group PL tag in the BAM header
    seq_platform: String = ""

    // Whether to ignore the GTF in STAR alignment
    star_ignore_sjdbgtf: Boolean

    // The maximum amount of RAM to use for sorting the BAM file in STAR. Should by in bits. Setting this value to `0` will use the default amount of STAR.
    star_limit_bam_sort_ram: Integer = 0
}

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden,
        params.no_cosmic,
        params.dfam_version,
        params.species,
        [
            'dfam_hmm': params.dfam_hmm,
            'dfam_h3f': params.dfam_h3f,
            'dfam_h3i': params.dfam_h3i,
            'dfam_h3m': params.dfam_h3m,
            'dfam_h3p': params.dfam_h3p,
        ],
        params.pfam_version,
        params.pfam_file,
        params.genomes,
        params.genome
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
        // channels
        PIPELINE_INITIALISATION.out.samplesheet,
        // lists
        tools,
        // strings
        params.seq_center,
        params.seq_platform,
        params.fasta,
        params.fai,
        params.gtf,
        params.genome_gencode_version,
        params.genome,
        params.hgnc_ref,
        params.hgnc_date,
        params.rrna_intervals,
        params.refflat,
        params.salmon_index,
        params.starindex_ref,
        params.arriba_ref_blacklist,
        params.arriba_ref_cytobands,
        params.arriba_ref_known_fusions,
        params.arriba_ref_protein_domains,
        params.fusioncatcher_ref,
        params.starfusion_ref,
        params.fusionreport_ref,
        params.fusion_annot_lib,
        params.species,
        params.cosmic_username,
        params.cosmic_passwd,
        // numbers
        params.min_trimmed_reads,
        params.trim_tail_fusioncatcher,
        params.tools_cutoff,
        // booleans
        params.references_only,
        params.save_trimmed_fail,
        params.save_merged,
        params.skip_qc,
        params.star_ignore_sjdbgtf,
        params.skip_vis,
        params.skip_vcf,
        params.no_cosmic,
        // paths
        params.adapter_fasta,
        params.arriba_fusions,
        params.starfusion_fusions,
        params.fusioncatcher_fusions,
        params.fusioninspector_fusions,
        params.whitelist,
        params.outdir,
        params.multiqc_methods_description,
        params.multiqc_config,
        params.multiqc_logo,
        params.pfam_file,
        params.dfam_hmm,
        params.dfam_h3f,
        params.dfam_h3i,
        params.dfam_h3m,
        params.dfam_h3p,
        params.annot_filter_url,
        params.ctatsplicing_cancer_introns,
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
        RNAFUSION.out.multiqc_report,
        params.max_multiqc_email_size
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
