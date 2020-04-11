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

    nextflow run nf-core/rnafusion/download-references.nf -profile [PROFILE] [OPTIONS] --outdir /path/to/output

    Mandatory arguments:
      --reference_release [int]     Release number of Ensembl reference for FASTA and GTF
                                    example: 97 -> ftp://ftp.ensembl.org/pub/release-97
      --outdir [path]               Output directory for downloading
      -profile [str]                Configuration profile [https://github.com/nf-core/configs]
      
    Options:
      --download_all [bool]         Download all references
      --arriba [bool]               Download Arriba references
      --star_fusion [bool]          Build STAR-Fusion references from FASTA ANF GTF
      --fusioncatcher [bool]        Download Fusioncatcher references
      --ericscript [bool]           Download Ericscript references 
      --fusion_report [bool]        Download databases for fusion-report
        --cosmic_usr [str]          [Required] COSMIC username
        --cosmic_passwd [str]       [Required] COSMIC password
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help) exit 0, helpMessage()

if (!params.reference_release) exit 1, "You did not specify the release number of reference!"
if (!params.outdir) exit 1, "Output directory not specified!"

params.running_tools = []
if (params.arriba || params.download_all) params.running_tools.add("Arriba")
if (params.star_fusion || params.download_all) params.running_tools.add("STAR-Fusion")
if (params.fusioncatcher || params.download_all) params.running_tools.add("Fusioncatcher")
if (params.ericscript || params.download_all) params.running_tools.add("Ericscript")
if (params.download_all) {
    params.running_tools.add('fusion-report')
    if (!params.cosmic_usr || !params.cosmic_passwd) exit 1, "Database credentials are required parameter!"
}

// Header log info
log.info nfcoreHeader()
def summary = [:]
summary['Pipeline Name']  = 'nf-core/rnafusion/download-references.nf'
summary['Pipeline Version'] = workflow.manifest.version
summary['References']       = params.running_tools.size() == 0 ? 'None' : params.running_tools.join(", ")
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

process download_base {
    publishDir "${params.outdir}/", mode: 'copy'
    
    output:
    file "Homo_sapiens.GRCh38_r${params.reference_release}.all.fa" into fasta
    file "Homo_sapiens.GRCh38_r${params.reference_release}.gtf" into gtf
    file "Homo_sapiens.GRCh38_${params.reference_release}.cdna.all.fa.gz" into transcript

    script:
    """
    wget ftp://ftp.ensembl.org/pub/release-${params.reference_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{1..22}.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-${params.reference_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{MT,X,Y}.fa.gz
    gunzip -c Homo_sapiens.GRCh38.dna.chromosome.* > Homo_sapiens.GRCh38_r${params.reference_release}.all.fa
    wget ftp://ftp.ensembl.org/pub/release-${params.reference_release}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${params.reference_release}.chr.gtf.gz -O Homo_sapiens.GRCh38_r${params.reference_release}.gtf.gz
    gunzip Homo_sapiens.GRCh38_r${params.reference_release}.gtf.gz
    wget ftp://ftp.ensembl.org/pub/release-${params.reference_release}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -O Homo_sapiens.GRCh38_${params.reference_release}.cdna.all.fa.gz
    """
}

process download_arriba {
    publishDir "${params.outdir}/arriba", mode: 'copy'
    
    when:
    params.arriba || params.download_all

    output:
    file '*'

    script:
    """
    wget -N https://github.com/suhrig/arriba/releases/download/v1.1.0/arriba_v1.1.0.tar.gz -O arriba_v1.1.0.tar.gz
    tar -xvzf arriba_v1.1.0.tar.gz && mv arriba_v1.1.0/database/* . && gunzip *.gz && rm -rf arriba_*
    """
}

process download_star_fusion {
    label 'process_high'
    publishDir "${params.outdir}/star-fusion", mode: 'copy'
    
    when:
    params.star_fusion || params.download_all

    input:
    file fasta
    file gtf

    output:
    file '*'

    script:
    """
    wget -N ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
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
        --annot_filter_rule /opt/conda/envs/star-fusion_v1.6.0/ctat-genome-lib-builder-2830cd708c5bb9353878ca98069427e83acdd625/AnnotFilterRuleLib/AnnotFilterRule.pm \\
        --fusion_annot_lib CTAT_HumanFusionLib.dat.gz \\
        --pfam_db Pfam-A.hmm \\
        --dfam_db homo_sapiens_dfam.hmm \\
        --CPU ${task.cpus}
    """
}

process download_fusioncatcher {
    publishDir "${params.outdir}/fusioncatcher", mode: 'copy'
    
    when:
    params.fusioncatcher || params.download_all

    output:
    file '*'

    script:
    """
    wget -N http://sourceforge.net/projects/fusioncatcher/files/data/human_v98.tar.gz.aa
    wget -N http://sourceforge.net/projects/fusioncatcher/files/data/human_v98.tar.gz.ab
    wget -N http://sourceforge.net/projects/fusioncatcher/files/data/human_v98.tar.gz.ac
    wget -N http://sourceforge.net/projects/fusioncatcher/files/data/human_v98.tar.gz.ad
    cat human_v98.tar.gz.* | tar xz
    rm human_v98.tar*
    """
}

process download_ericscript {
    publishDir "${params.outdir}/ericscript", mode: 'copy'
    
    when:
    params.ericscript || params.download_all

    output:
    file '*'

    script:
    """
    wget -N https://raw.githubusercontent.com/circulosmeos/gdown.pl/dfd6dc910a38a42d550397bb5c2335be2c4bcf54/gdown.pl
    chmod +x gdown.pl
    ./gdown.pl "https://drive.google.com/uc?export=download&confirm=qgOc&id=0B9s__vuJPvIiUGt1SnFMZFg4TlE" ericscript_db_homosapiens_ensembl84.tar.bz2
    tar jxf ericscript_db_homosapiens_ensembl84.tar.bz2
    rm gdown.pl ericscript_db_homosapiens_ensembl84.tar.bz2
    """
}

process download_databases {
    publishDir "${params.outdir}/databases", mode: 'copy'

    output:
    file '*'

    script:
    """
    fusion_report download --cosmic_usr "${params.cosmic_usr}" --cosmic_passwd "${params.cosmic_passwd}" .
    """
}

/*
 * Completion
 */
workflow.onComplete {
    log.info "[nf-core/rnafusion/download-references.nf] Pipeline Complete"
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