process FUSIONINSPECTOR_DEV {
    tag "$meta.id"
    label 'process_high'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "docker.io/trinityctat/fusioninspector:2.8.0-dev1"
    } else {
        container "docker.io/trinityctat/fusioninspector:2.8.0-dev1"
    }

    input:
    tuple val(meta), path(reads), path(fusion_list)
    path reference

    output:
    path "*"                , emit: output
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta = meta.single_end ? "--left_fq ${reads[0]}" : "--left_fq ${reads[0]} --right_fq ${reads[1]}"
    def args = task.ext.args ?: ''
    """
    FusionInspector \\
        --fusions $fusion_list \\
        --genome_lib ${reference} \\
        $fasta \\
        --CPU ${task.cpus} \\
        -O . \\
        --out_prefix $prefix \\
        --vis $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR-Fusion: \$(echo STAR-Fusion --version 2>&1 | grep -i 'version' | sed 's/STAR-Fusion version: //')
    END_VERSIONS
    """

    stub:
    """
    touch versions.yml
    touch FusionInspector.log
    """
}
