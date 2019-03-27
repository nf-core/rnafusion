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
      --star_fusion                 Download STAR-Fusion references [NCBI version by default]
      --star_fusion_ensembl         Download STAR-Fusion references [Ensebml, have to build manually]
      --fusioncatcher               Download Fusioncatcher references
      --ericscript                  Download Ericscript references 
      --pizzly                      Download pizzly references
      --databases                   Download databases for fusion-report
        --cosmic_usr                [Required] COSMIC username
        --cosmic_passwd             [Required] COSMIC password
      --igenomesIgnore              Download iGenome Homo Sapiens version NCBI/GRCh38.
                                    Ignored on default
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
if (params.download_all) {
    params.running_tools.add("All")
}
if (params.igenome) {
    params.running_tools.add("iGenome")
}
if (params.star_fusion) {
    params.running_tools.add("STAR-Fusion NCBI version")
}
if (params.star_fusion_ensembl) {
    params.running_tools.add("STAR-Fusion Ensembl version")
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
if (params.databases) {
    if (!params.cosmic_usr && params.cosmic_passwd) {
        exit 1, "Database credentials are required parameter!"
    }
    params.running_tools.add('Databases')
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
    params.star_fusion || params.download_all

    output:
    file '*'

    script:
    """
    wget -N https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz -O GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz
    tar -xvzf GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz && rm GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz
    """
}

process download_star_fusion_ensembl {
    publishDir "${params.outdir}/star_fusion_ensembl_ref", mode: 'copy'
    
    when:
    params.star_fusion_ensembl || params.download_all

    output:
    file '*'

    script:
    """
    wget ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{1..22}.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{MT,X,Y}.fa.gz
    gunzip -c Homo_sapiens.GRCh38.dna.chromosome.* > Homo_sapiens.GRCh38_r77.all.fa
    wget ftp://ftp.ensembl.org/pub/release-77/gtf/homo_sapiens/Homo_sapiens.GRCh38.77.chr.gtf.gz
    wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
    gunzip Pfam-A.hmm.gz && hmmpress Pfam-A.hmm
    prep_genome_lib.pl \\
        --genome_fa Homo_sapiens.GRCh38_r77.all.fa \\
        --gtf Homo_sapiens.GRCh38.77.chr.gtf \\
        --pfam_db Pfam-A.hmm \\
        --CPU 10
    """
}

process download_fusioncatcher {
    publishDir "${params.outdir}/fusioncatcher_ref", mode: 'copy'
    
    when:
    params.fusioncatcher || params.download_all

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

process download_pizzly {
    publishDir "${params.outdir}/pizzly_ref", mode: 'copy'
    
    when:
    params.pizzly || params.download_all

    output:
    file '*'

    script:
    """
    wget -N ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
    wget -N ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz && gunzip Homo_sapiens.GRCh38.94.gtf.gz
    """
}

process download_databases {
    publishDir "${params.outdir}/databases", mode: 'copy'
    
    when:
    params.databases || params.download_all

    output:
    file '*'

    script:
    """
    fusion_report download --cosmic_usr ${params.cosmic_usr} --cosmic_passwd ${cosmic_passwd} .
    """
}

process download_igenome {
    publishDir "${params.outdir}/igenome", mode: 'copy'
    
    when:
    !params.igenomesIgnore

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