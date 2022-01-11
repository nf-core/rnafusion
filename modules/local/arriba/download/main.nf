process ARRIBA_DOWNLOAD {
    tag "arriba"
    label 'process_low'
    publishDir "${params.genomes_base}",
        mode: params.publishDir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }

    conda (params.enable_conda ? "bioconda::gnu-wget=1.18" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h5bf99c6_5"
    } else {
        container "quay.io/biocontainers/gnu-wget:1.18--h5bf99c6_5"
    }

    output:
    path "versions.yml"            , emit: versions
    path "arriba_v2.1.0/database/*" , emit: reference

    script:
    """
    wget https://github.com/suhrig/arriba/releases/download/v2.1.0/arriba_v2.1.0.tar.gz -O arriba_v2.1.0.tar.gz
    tar -xzvf arriba_v2.1.0.tar.gz
    rm arriba_v2.1.0.tar.gz

    echo \$(wget -V 2>&1) | grep "GNU Wget" | cut -d" " -f3 > versions.yml
    """
}
