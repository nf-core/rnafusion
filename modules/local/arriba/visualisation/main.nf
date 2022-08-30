process ARRIBA_VISUALISATION {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::arriba=2.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/arriba:2.1.0--h3198e80_1' :
        'quay.io/biocontainers/arriba:2.1.0--h3198e80_1' }"

    input:
    tuple val(meta), path(bam), path(bai), path(fusions)
    path reference
    path gtf

    output:
    tuple val(meta), path("*.pdf")          , emit: pdf
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    draw_fusions.R \\
        --fusions=$fusions \\
        --alignments=$bam \\
        --output=${prefix}.pdf \\
        --annotation=${gtf} \\
        --cytobands=${reference}/${args} \\
        --proteinDomains=${reference}/${args2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arriba: \$(arriba -h | grep 'Version:' 2>&1 |  sed 's/Version:\s//')
    END_VERSIONS
    """
}
