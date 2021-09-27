// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PIZZLY_DOWNLOAD {
    tag "pizzly"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::kallisto=0.46.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/kallisto:0.46.2--h4f7b962_1"
    } else {
        container "quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1"
    }

    input:
    path transcript

    output:
    path "*.version.txt", emit: version
    path "index.idx"    , emit: reference

    script:
    def software = getSoftwareName(task.process)
    """
    kallisto index \\
        -i index.idx \\
        $options.args \\
        $transcript

    echo \$(kallisto 2>&1) | sed 's/^kallisto //; s/Usage.*\$//' > ${software}.version.txt
    """
}