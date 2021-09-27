// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARRIBA_DOWNLOAD {
    tag "arriba"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::gnu-wget=1.18" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h5bf99c6_5"
    } else {
        container "quay.io/biocontainers/gnu-wget:1.18--h5bf99c6_5"
    }

    output:
    path "*.version.txt"           , emit: version
    path "arriba_v2.1.0/database/*", emit: reference

    script:
    def software = getSoftwareName(task.process)
    """
    wget https://github.com/suhrig/arriba/releases/download/v2.1.0/arriba_v2.1.0.tar.gz -O arriba_v2.1.0.tar.gz
    tar -xzvf arriba_v2.1.0.tar.gz
    rm arriba_v2.1.0.tar.gz

    echo \$(wget -V 2>&1) | grep "GNU Wget" | cut -d" " -f3 > ${software}.version.txt
    """
}