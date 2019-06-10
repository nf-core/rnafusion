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

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    
    Usage:

    The typical command for downloading singularity images is as follows:

    nextflow run nf-core/rnafusion/download-singularity-img.nf -profile [PROFILE] [OPTIONS] --outdir /path/to/output

    By default main image is downloaded.

    Mandatory arguments:
      --outdir                      Output directory for downloading
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: standard, conda, docker, singularity, awsbatch, test
      
    Options:
      --download_all                Download all images
      --arriba                      Download Arriba image
      --ericscript                  Download Ericscript image 
      --fusioncatcher               Download Fusioncatcher image
      --fusion_inspector            Download Fusion-Inspector image
      --pizzly                      Download Pizzly image
      --squid                       Download Squid image
      --star_fusion                 Download STAR-Fusion image
      
    """.stripIndent()
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
if (params.download_all) {
    params.running_tools.add("All")
}
if (params.arriba) {
    params.running_tools.add("Arriba")
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
log.info nfcoreHeader()
def summary = [:]
summary['Pipeline Name']  = 'nf-core/rnafusion/download-singularity-img.nf'
summary['Pipeline Version'] = workflow.manifest.version
summary['Tool images']      = params.running_tools.size() == 0 ? 'None' : params.running_tools.join(", ")
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Launch dir']   = workflow.launchDir
summary['Working dir']  = workflow.workDir
summary['Script dir']   = workflow.projectDir
summary['User']         = workflow.userName
summary['Config Profile'] = workflow.profile
if(params.config_profile_description) summary['Config Description'] = params.config_profile_description
if(params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if(params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if(workflow.profile == 'awsbatch'){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "\033[2m----------------------------------------------------\033[0m"

process download_base_image {
    publishDir "${params.outdir}", mode: 'copy'

    when:
    params.download_all

    output:
    file "rnafusion_v${workflow.manifest.version}.img"

    script:
    """
    singularity pull --name "rnafusion_v${workflow.manifest.version}.img" docker://nfcore/rnafusion:${workflow.manifest.version}
    """
}

process download_arriba {
    publishDir "${params.outdir}", mode: 'copy'
    
    when:
    params.arriba || params.download_all

    output:
    file "rnafusion_arriba_v${params.arriba_version}.img"

    script:
    """
    singularity pull --name "rnafusion_arriba_v${params.arriba_version}.img" docker://nfcore/rnafusion:arriba_v${params.arriba_version}
    """
}

process download_ericscript {
    publishDir "${params.outdir}", mode: 'copy'
    
    when:
    params.ericscript || params.download_all

    output:
    file "rnafusion_ericscript_v${params.ericscript_version}.img"

    script:
    """
    singularity pull --name "rnafusion_ericscript_v${params.ericscript_version}.img" docker://nfcore/rnafusion:ericscript_v${params.ericscript_version}
    """
}

process download_fusioncatcher {
    publishDir "${params.outdir}", mode: 'copy'
    
    when:
    params.fusioncatcher || params.download_all

    output:
    file "rnafusion_fusioncatcher_v${params.fusioncatcher_version}.img"

    script:
    """
    singularity pull --name "rnafusion_fusioncatcher_v${params.fusioncatcher_version}.img" docker://nfcore/rnafusion:fusioncatcher_v${params.fusioncatcher_version}
    """
}

process download_fusion_inspector {
    publishDir "${params.outdir}", mode: 'copy'
    
    when:
    params.fusion_inspector || params.download_all

    output:
    file "rnafusion_fusion-inspector_v${params.fusion_inspector_version}.img"

    script:
    """
    singularity pull --name "rnafusion_fusion-inspector_v${params.fusion_inspector_version}.img" docker://nfcore/rnafusion:fusion-inspector_v${params.fusion_inspector_version}
    """
}

process download_pizzly {
    publishDir "${params.outdir}", mode: 'copy'
    
    when:
    params.pizzly || params.download_all

    output:
    file "rnafusion_pizzly_v${params.pizzly_version}.img"

    script:
    """
    singularity pull --name "rnafusion_pizzly_v${params.pizzly_version}.img" docker://nfcore/rnafusion:pizzly_v${params.pizzly_version}
    """
}

process download_squid {
    publishDir "${params.outdir}", mode: 'copy'
    
    when:
    params.squid || params.download_all

    output:
    file "rnafusion_squid_v${params.squid_version}.img"

    script:
    """
    singularity pull --name "rnafusion_squid_v${params.squid_version}.img" docker://nfcore/rnafusion:squid_v${params.squid_version}
    """
}

process download_star_fusion {
    publishDir "${params.outdir}", mode: 'copy'
    
    when:
    params.star_fusion || params.download_all

    output:
    file "rnafusion_star-fusion_v${params.star_fusion_version}.img"

    script:
    """
    singularity pull --name "rnafusion_star-fusion_v${params.star_fusion_version}.img" docker://nfcore/rnafusion:star-fusion_v${params.star_fusion_version}
    """
}

/*
 * Completion
 */
workflow.onComplete {
    log.info "[nf-core/rnafusion] Pipeline Complete"
}

def nfcoreHeader(){
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """    ${c_dim}----------------------------------------------------${c_reset}
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/rnafusion v${workflow.manifest.version}${c_reset}
    ${c_dim}----------------------------------------------------${c_reset}
    """.stripIndent()
}