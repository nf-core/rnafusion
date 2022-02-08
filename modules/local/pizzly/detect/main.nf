process PIZZLY {
    tag "pizzly"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::kallisto=0.46.2 bioconda::pizzly==0.37.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/kallisto:0.46.2--h4f7b962_1"
    } else {
        container "quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1"
    }

    input:
    tuple val(meta), path(txt)
    transcript
    gtf

    output:
    path "versions.yml"                                     , emit: versions
    tuple val(meta), path "${prefix}.pizzly.txt"            , emit: fusions
    tuple val(meta), path "${prefix}.pizzly.unfiltered.json", emit: fusions_unfiltered

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pizzly \\
        $args \\
        --gtf $gtf \\
        --fasta $transcript \\
        --output ${prefix}.pizzly output/fusion.txt

    pizzly_flatten_json.py ${prefix}.pizzly.json ${prefix}.pizzly.txt

    echo \$(kallisto 2>&1) | sed 's/^kallisto //; s/Usage.*\$//' > versions.yml
    """
}
