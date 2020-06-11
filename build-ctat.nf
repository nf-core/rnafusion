#!/usr/bin/env nextflow
/*
================================================================================
                                nf-core/rnafusion
================================================================================
nf-core/rnafusion:
 RNA-seq analysis pipeline for detection gene-fusions
--------------------------------------------------------------------------------
 @Homepage
 https://nf-co.re/rnafusion
--------------------------------------------------------------------------------
 @Documentation
 https://nf-co.re/rnafusion/docs
--------------------------------------------------------------------------------
 @Repository
 https://github.com/nf-core/rnafusion
--------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Usage:

    The typical command for downloading references is as follows:

    nextflow run nf-core/rnafusion/build-ctat.nf -profile [PROFILE] [OPTIONS] --outdir /path/to/output

    Mandatory arguments:
      --fasta [file]                Path to fasta reference
      --gtf [file]                  Path to GTF annotation
      --outdir [path]               Output directory for downloading
      -profile [str]                Configuration profile [https://github.com/nf-core/configs]
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help) exit 0, helpMessage()

params.fasta = params.genome ? params.genomes[params.genome].fasta ?: null : null
params.gtf = params.genome ? params.genomes[params.genome].gtf ?: null : null

ch_fasta = Channel.value(file(params.fasta)).ifEmpty{exit 1, "Fasta file not found: ${params.fasta}"}
ch_gtf = Channel.value(file(params.gtf)).ifEmpty{exit 1, "GTF annotation file not found: ${params.gtf}"}

if (!params.outdir) exit 1, "Output directory not specified!"

// Header log info
log.info nfcoreHeader()
def summary = [:]
summary['Pipeline Name']  = 'nf-core/rnafusion/build-ctat.nf'
summary['Pipeline Version'] = workflow.manifest.version
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
summary['Output dir']   = params.outdir
summary['User']         = workflow.userName
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "\033[2m----------------------------------------------------\033[0m"

// Check the hostnames against configured profiles
checkHostname()

/*
================================================================================
                                  DOWNLOAD
================================================================================
*/

process star_fusion {
    label 'process_high'
    label 'process_long'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file(fasta) from ch_fasta
    file(gtf) from ch_gtf

    output:
    file '*'

    script:
    """
    wget -N ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_releases/Pfam-A.hmm.gz
    gunzip Pfam-A.hmm.gz && hmmpress Pfam-A.hmm

    wget https://github.com/FusionAnnotator/CTAT_HumanFusionLib/releases/download/v0.2.0/fusion_lib.Mar2019.dat.gz -O CTAT_HumanFusionLib.dat.gz
    
    # Dfam
    wget https://www.dfam.org/releases/Dfam_3.1/infrastructure/dfamscan/homo_sapiens_dfam.hmm
    wget https://www.dfam.org/releases/Dfam_3.1/infrastructure/dfamscan/homo_sapiens_dfam.hmm.h3f
    wget https://www.dfam.org/releases/Dfam_3.1/infrastructure/dfamscan/homo_sapiens_dfam.hmm.h3i
    wget https://www.dfam.org/releases/Dfam_3.1/infrastructure/dfamscan/homo_sapiens_dfam.hmm.h3m
    wget https://www.dfam.org/releases/Dfam_3.1/infrastructure/dfamscan/homo_sapiens_dfam.hmm.h3p

    export TMPDIR=/tmp
    prep_genome_lib.pl \\
        --genome_fa ${fasta} \\
        --gtf ${gtf} \\
        --annot_filter_rule /opt/conda/envs/nf-core-rnafusion-star-fusion_1.8.1/ctat-genome-lib-builder-2830cd708c5bb9353878ca98069427e83acdd625/AnnotFilterRuleLib/AnnotFilterRule.pm \\
        --fusion_annot_lib CTAT_HumanFusionLib.dat.gz \\
        --pfam_db Pfam-A.hmm \\
        --dfam_db homo_sapiens_dfam.hmm \\
        --CPU ${task.cpus}
    """
}

/*
 * Completion
 */
workflow.onComplete {
    log.info "[nf-core/rnafusion/build-ctat.nf] Pipeline Complete"
}

def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/rnafusion v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
