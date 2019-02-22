#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/rnafusion
========================================================================================
 nf-core/rnafusion Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/rnafusion
----------------------------------------------------------------------------------------
*/

nfcore_logo = """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

nf-core/rnafusion v${workflow.manifest.version}
======================================================="""

def helpMessage() {
    nfcore_help = """
    Usage:

    The typical command for downloading references is as follows:

    nextflow run nf-core/rnafusion download-references.nf [OPTIONS] --outdir /path/to/output

    Mandatory arguments:
      --outdir                      Output directory for downloading
      
    Options:
      --all                         Download all references except iGenome
      --star_fusion                 Download STAR-Fusion references
      --fusioncatcher               Download Fusioncatcher references
      --ericscript                  Download Ericscript references 
      --pizzly                      Download pizzly references
      --igenome                     Download iGenome Homo Sapiens version NCBI/GRCh38
    """.stripIndent()
    log.info "${nfcore_logo}${nfcore_help}"
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

params.running_tools = []
if (!params.outdir) {
    exit 1, "Output directory not specified!"
}
if (params.all) {
    params.running_tools.add("All")
}
if (params.igenome) {
    params.running_tools.add("iGenome")
}
if (params.star_fusion) {
    params.running_tools.add("STAR-Fusion")
}
if (params.fusioncatcher) {
    params.running_tools.add("Fusioncatcher")
}
if (params.ericscript) {
    params.running_tools.add("Ericscript")
}
if (params.pizzly) {
    params.running_tools.add("Pizzly")
}

// Header log info
log.info nfcore_logo
def summary = [:]
summary['Pipeline Name']  = 'nf-core/rnafusion/download-references.nf'
summary['Pipeline Version'] = workflow.manifest.version
summary['References']       = params.running_tools.size() == 0 ? 'None' : params.running_tools.join(", ")
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

process download_star_fusion {
    publishDir "${params.outdir}/star_fusion_ref", mode: 'copy'
    
    when:
    params.star_fusion || params.all

    output:
    file '*'

    script:
    """
    wget -N https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz -O GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz
    tar -xvzf GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz && rm GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz
    """
}

process download_fusioncatcher {
    publishDir "${params.outdir}/fusioncatcher_ref", mode: 'copy'
    
    when:
    params.fusioncatcher || params.all

    output:
    file '*'

    script:
    """
    wget -N http://sourceforge.net/projects/fusioncatcher/files/data/human_v90.tar.gz.aa
    wget -N http://sourceforge.net/projects/fusioncatcher/files/data/human_v90.tar.gz.ab
    wget -N http://sourceforge.net/projects/fusioncatcher/files/data/human_v90.tar.gz.ac
    wget -N http://sourceforge.net/projects/fusioncatcher/files/data/human_v90.tar.gz.ad
    cat human_v90.tar.gz.* | tar xz
    rm human_v90.tar*
    """
}

process download_ericscript {
    publishDir "${params.outdir}/ericscript_ref", mode: 'copy'
    
    when:
    params.ericscript || params.all

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

process download_pizzly {
    publishDir "${params.outdir}/pizzly_ref", mode: 'copy'
    
    when:
    params.pizzly || params.all

    output:
    file '*'

    script:
    """
    wget -N ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
    wget -N ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz && gunzip Homo_sapiens.GRCh38.94.gtf.gz
    """
}

process download_igenome {
    publishDir "${params.outdir}/igenome", mode: 'copy'
    
    when:
    params.igenome

    output:
    file '*'

    script:
    """
    aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/NCBI/GRCh38/ .
    """
}

/*
 * Completion
 */
workflow.onComplete {
    log.info "[nf-core/rnafusion] Pipeline Complete"
}