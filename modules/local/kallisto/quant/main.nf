process KALLISTO_QUANT {
    tag "kalliso_quant"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::kallisto=0.46.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kallisto:0.46.2--h4f7b962_1' :
        'quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1' }"


    input:
    tuple val(meta), path(reads)
    path index

    output:
    path "versions.yml"                        , emit: versions
    tuple val(meta), path("fusion.txt")        , emit: txt

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

    touch versions.yml

    """
}
    // cat <<-END_VERSIONS > versions.yml
    // "${task.process}":
    //     kallisto: \$(kallisto | sed 's/^kallisto //; s/Usage.*\$//')
    // END_VERSIONS
    // echo \$(kallisto 2>&1) | sed 's/^kallisto //; s/Usage.*\$//' > versions.yml
