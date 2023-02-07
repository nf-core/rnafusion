process KALLISTO_QUANT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::kallisto=0.46.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kallisto:0.46.2--h4f7b962_1' :
        'quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1' }"


    input:
    tuple val(meta), path(reads)
    path index

    output:
    path "versions.yml"                        , emit: versions
    tuple val(meta), path("*fusions.txt")      , emit: txt

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    kallisto quant \
        -t $task.cpus \
        -i $index \
        --fusion \
        -o . \
        $reads
    mv fusion.txt ${prefix}.kallisto_quant.fusions.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$(echo \$(kallisto 2>&1) | sed 's/^kallisto //; s/Usage.*\$//')
    END_VERSIONS
    """
}

