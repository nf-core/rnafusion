// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FUSIONCATCHER_DOWNLOADGENOME {
    tag '$organism'
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'resources', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::fusioncatcher=1.33" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/fusioncatcher:1.33--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/fusioncatcher:1.33--hdfd78af_1"
    }

    input:
    val organism

    output:
    path "fusioncatcher-genome"     , emit: reference
    path "*.version.txt"            , emit: version

    script:
    def software = getSoftwareName(task.process)

    """
    fusioncatcher-build \\
        -g $organism \\
        -o fusioncatcher-genome \\
        $options.args

    fusioncatcher --version | sed 's/fusioncatcher.py //' > ${software}.version.txt
    """
}
