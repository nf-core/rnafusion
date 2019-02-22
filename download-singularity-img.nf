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

    The typical command for downloading singularity images is as follows:

    nextflow run nf-core/rnafusion download-singularity-img.nf [OPTIONS] --outdir /path/to/output

    By default main image is downloaded.

    Mandatory arguments:
      --outdir                      Output directory for downloading
      
    Options:
      --all                         Download all images
      --star_fusion                 Download STAR-Fusion image
      --fusioncatcher               Download Fusioncatcher image
      --ericscript                  Download Ericscript image 
      --pizzly                      Download Pizzly image
      --squid                       Download Squid image
      --fusion_inspector            Download Fusion-Inspector image
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

params.running_tools = ["nf-core/rnafusion/${workflow.manifest.version}"]
if (!params.outdir) {
    exit 1, "Output directory not specified!"
}
if (params.all) {
    params.running_tools.add("All")
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
if (params.squid) {
    params.running_tools.add("Squid")
}
if (params.fusion_inspector) {
    params.running_tools.add("Fusion-Inspector")
}

// Header log info
log.info nfcore_logo
def summary = [:]
summary['Pipeline Name']  = 'nf-core/rnafusion/download-singularity-img.nf'
summary['Pipeline Version'] = workflow.manifest.version
summary['Tool images']      = params.running_tools.size() == 0 ? 'None' : params.running_tools.join(", ")
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

process download_base_image {
    publishDir "${params.outdir}", mode: 'copy'

    output:
    file '*'

    script:
    """
    singularity pull --name rnafusion_v${workflow.manifest.version}.img docker://nfcore/rnafusion:${workflow.manifest.version}
    """
}

process download_star_fusion {
    publishDir "${params.outdir}", mode: 'copy'
    
    when:
    params.star_fusion || params.all

    output:
    file '*'

    script:
    """
    singularity pull --name rnafusion_star-fusion_v${workflow.manifest.version}.img docker://nfcore/rnafusion:star-fusion_v${workflow.manifest.version}
    """
}

process download_fusioncatcher {
    publishDir "${params.outdir}", mode: 'copy'
    
    when:
    params.fusioncatcher || params.all

    output:
    file '*'

    script:
    """
    singularity pull --name rnafusion_fusioncatcher_v${workflow.manifest.version}.img docker://nfcore/rnafusion:fusioncatcher_v${workflow.manifest.version}
    """
}

process download_ericscript {
    publishDir "${params.outdir}", mode: 'copy'
    
    when:
    params.ericscript || params.all

    output:
    file '*'

    script:
    """
    singularity pull --name rnafusion_ericscript_v${workflow.manifest.version}.img docker://nfcore/rnafusion:ericscript_v${workflow.manifest.version}
    """
}

process download_pizzly {
    publishDir "${params.outdir}", mode: 'copy'
    
    when:
    params.pizzly || params.all

    output:
    file '*'

    script:
    """
    singularity pull --name rnafusion_pizzly_v${workflow.manifest.version}.img docker://nfcore/rnafusion:pizzly_v${workflow.manifest.version}
    """
}

process download_squid {
    publishDir "${params.outdir}", mode: 'copy'
    
    when:
    params.squid || params.all

    output:
    file '*'

    script:
    """
    singularity pull --name rnafusion_squid_v${workflow.manifest.version}.img docker://nfcore/rnafusion:squid_v${workflow.manifest.version}
    """
}

process download_fusion_inspector {
    publishDir "${params.outdir}", mode: 'copy'
    
    when:
    params.fusion_inspector || params.all

    output:
    file '*'

    script:
    """
    singularity pull --name rnafusion_fusion-inspector_v${workflow.manifest.version}.img docker://nfcore/rnafusion:fusion-inspector_v${workflow.manifest.version}
    """
}

/*
 * Completion
 */
workflow.onComplete {
    log.info "[nf-core/rnafusion] Pipeline Complete"
}