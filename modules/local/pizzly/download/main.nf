process PIZZLY_DOWNLOAD {
    tag "pizzly"
    label 'process_medium'
    publishDir "${params.genomes_base}",
        mode: params.publishDir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }

    conda (params.enable_conda ? "bioconda::kallisto=0.46.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/kallisto:0.46.2--h4f7b962_1"
    } else {
        container "quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1"
    }

    input:
    path transcript

    output:
    path "versions.yml" , emit: versions
    path "index.idx"    , emit: reference

    script:
    """
    kallisto index \\
        -i index.idx \\
        $options.args \\
        $transcript

    echo \$(kallisto 2>&1) | sed 's/^kallisto //; s/Usage.*\$//' > versions.yml
    """
}
