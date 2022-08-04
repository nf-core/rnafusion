process REFORMAT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bbmap=38.90" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:38.90--he522d1c_1' :
        'quay.io/biocontainers/bbmap:38.90--he522d1c_1' }"


    input:
    tuple val(meta), path(reads)

    output:
    path "versions.yml"                                             , emit: versions
    tuple val(meta), path("*trimmed.fq.gz")                         , emit: reads_out


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def in1 = "in=${reads[0]}"
    def in2 = meta.single_end ? "" : "in=${reads[1]}"
    def out1 ="out=${prefix}_R1_trimmed.fq.gz"
    def out2 =meta.single_end ? "" : "out=${prefix}_R2_trimmed.fq.gz"

    """
    reformat.sh $in1 $out1 $args
    reformat.sh $in2 $out2 $args2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        reformat.sh: \$(echo \$(reformat.sh --version 2>&1)| sed -e "s/BBMap version //g" )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out1 ="out=${prefix}_R1_trimmed.fq.gz"
    def out2 =meta.single_end ? "" : "out=${prefix}_R2_trimmed.fq.gz"
    """
    touch $out1
    touch $out2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        reformat.sh: \$(echo \$(reformat.sh --version 2>&1)| sed -e "s/BBMap version //g" )
    END_VERSIONS
    """
}
