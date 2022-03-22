process PIZZLY {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::kallisto=0.46.2 bioconda::pizzly==0.37.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pizzly:0.37.3--py36_2"
    } else {
        container "quay.io/biocontainers/pizzly:0.37.3--h470a237_3"
    }

    input:
    tuple val(meta), path(txt)
    path transcript
    path gtf

    output:
    path "versions.yml"                              , emit: versions
    tuple val(meta), path('*pizzly.txt')            , emit: fusions
    tuple val(meta), path('*unfiltered.json')        , emit: fusions_unfiltered

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pizzly \\
        $args \\
        --gtf $gtf \\
        --fasta $transcript \\
        --output ${prefix}.pizzly $txt

    pizzly_flatten_json.py ${prefix}.pizzly.json ${prefix}.pizzly.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pizzly: \$(pizzly --version | grep pizzly | sed -e "s/pizzly version: //g")
    END_VERSIONS
    """
}
