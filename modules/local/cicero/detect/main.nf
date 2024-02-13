process CICERO {
    tag "$meta.id"
    label 'process_medium'

    container "docker.io/rannickscilifelab/cicero:1.9.6"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta), path(junctions)
    tuple val(meta), path(refdir)

    output:
    tuple val(meta), path("./*/*.final_fusions.txt") , emit: fusions
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    Cicero.sh \\
        -n $task.cpus \\
        -b $bam \\
        -r $refdir \\
        -j $junctions \\
        -o . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cicero: 1.9.6
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}/bla/${prefix}.final_fusions.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cicero: 1.9.6
    END_VERSIONS
    """
}
