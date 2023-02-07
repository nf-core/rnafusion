process PIZZLY_DOWNLOAD {
    tag "pizzly"
    label 'process_medium'

    conda "bioconda::kallisto=0.46.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kallisto:0.46.2--h4f7b962_1' :
        'quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1' }"

    input:
    path transcript

    output:
    path "versions.yml" , emit: versions
    path "index.idx"    , emit: reference

    script:
    def args = task.ext.args ?: ''
    """
    kallisto index \\
        -i index.idx \\
        $args \\
        $transcript

    echo \$(kallisto 2>&1) | sed 's/^kallisto //; s/Usage.*\$//' > versions.yml
    """

    stub:
    """
    touch index.idx
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$(echo \$(kallisto 2>&1) | sed 's/^kallisto //; s/Usage.*\$//')
    END_VERSIONS
    """
}
