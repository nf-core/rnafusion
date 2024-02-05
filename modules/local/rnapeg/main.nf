process RNAPEG {
    tag "$meta.id"
    label 'process_medium'

    container "docker.io/rannickscilifelab/rnapeg:2.7.7"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(refflat)

    output:
    tuple val(meta), path("*/junctions.txt")         , emit: junctions
    path "versions.yml"                              , emit: versions

    // when:
    // task.ext.when == null || task.ext.when

    // script:
    // def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    // """
    // RNApeg.sh \\
    //     -g GRCh38_no_alt \\
    //     -b $bam \\
    //     -f $fasta \\
    //     -r $refflat \\
    //     $args

    // cat <<-END_VERSIONS > versions.yml
    // "${task.process}":
    //     rnapeg: 2.7.7
    // END_VERSIONS
    // """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p results/${prefix}
    touch results/${prefix}/junctions.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnapeg: 2.7.7
    END_VERSIONS
    """
}
