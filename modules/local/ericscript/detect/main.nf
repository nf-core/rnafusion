process ERICSCRIPT {
    tag "eriscript"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::ericscript=0.5.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ericscript:0.5.5--pl5.22.0r3.4.1_1"
    } else {
        container "quay.io/biocontainers/ericscript:0.5.5--hdfd78af_5"
    }

    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path("*.results.filtered.tsv"), emit: fusions
    tuple val(meta), path("*.results.total.tsv")   , emit: fusions_total
    path "versions.yml"                           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    ericscript.pl \\
        -db $reference \\
        -name ${prefix} \\
        -p $task.cpus \\
        -o . \\
        $reads \\
        $args

    echo \$(wget -V 2>&1) | grep "GNU Wget" | cut -d" " -f3 > versions.yml

    """

    stub:
    """
    touch ${prefix}.results.filtered.tsv
    touch ${prefix}.results.total.tsv
    touch versions.yml
    """
}
