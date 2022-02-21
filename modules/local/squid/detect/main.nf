
process SQUID {
    tag "squid"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::squid=1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/nfcore/rnafusion:squid_1.5-star2.7.1a' :
        'docker.io/nfcore/rnafusion:squid_1.5-star2.7.1a' }"



    input:
    tuple val(meta), path(bam), path(chimeric_bam)

    output:
    tuple val(meta), path("*_sv.txt"), emit: fusions
    path  "versions.yml"          , emit: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    squid -b $bam -c $chimeric_bam -o ${prefix}_fusions

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        squid: \$(echo \$(squid --version 2>&1) | sed 's/v//')
    END_VERSIONS
    """
}
