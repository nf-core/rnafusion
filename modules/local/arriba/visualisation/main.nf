process ARRIBA_VISUALISATION {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::arriba=2.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/arriba:2.3.0--haa8aa89_0' :
        'quay.io/biocontainers/arriba:2.3.0--haa8aa89_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(fusions)
    path reference
    path gtf
    path protein_domains
    path cytobands

    output:
    tuple val(meta), path("*.pdf")          , emit: pdf
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def cytobands = cytobands ? " --cytobands=$cytobands" : ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def protein_domains = protein_domains ? "-p $protein_domains" : ""
    """
    draw_fusions.R \\
        --fusions=$fusions \\
        --alignments=$bam \\
        --output=${prefix}.pdf \\
        --annotation=${gtf} \\
        $cytobands \\
        $protein_domains \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arriba: \$(arriba -h | grep 'Version:' 2>&1 |  sed 's/Version:\s//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pdf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arriba: \$(arriba -h | grep 'Version:' 2>&1 |  sed 's/Version:\s//')
    END_VERSIONS
    """
}
