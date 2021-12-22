process SAMTOOLS_FAIDX {
    tag "$fasta"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publishDir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publishDir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? 'bioconda::samtools=1.13' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.13--h8c37831_0"
    } else {
        container "quay.io/biocontainers/samtools:1.13--h8c37831_0"
    }

    input:
    path fasta

    output:
    path "*.fai"        , emit: fai
    path "versions.yml" , emit: versions

    script:
    """
    samtools faidx $fasta
    echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' > versions.yml
    """
}
