
process SQUID {
    tag "squid"
    label 'process_medium'

    conda "bioconda::squid=1.5"
    container "docker.io/nfcore/rnafusion:squid_1.5-star2.7.1a"



    input:
    tuple val(meta), path(bam), path(chimeric_bam)

    output:
    tuple val(meta), path("*sv.txt") , emit: fusions
    path  "versions.yml"             , emit: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    squid -b $bam -c $chimeric_bam -o ${prefix}.squid.fusions

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        squid: \$(echo \$(squid --version 2>&1) | sed 's/v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.squid.fusions_sv.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        squid: \$(echo \$(squid --version 2>&1) | sed 's/v//')
    END_VERSIONS
    """
}
