process CONVERT2BED {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bedops=2.4.41"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedops:2.4.41--h9f5acd7_0' :
        'quay.io/biocontainers/bedops:2.4.41--h9f5acd7_0' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.bed") , emit: bed

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    convert2bed -i gtf < $gtf > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        convert2bed: \$(convert2bed --version | grep vers | sed 's/^.*.version:  //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        convert2bed: \$(convert2bed --version | grep vers | sed 's/^.*.version:  //')
    END_VERSIONS
    """
}
