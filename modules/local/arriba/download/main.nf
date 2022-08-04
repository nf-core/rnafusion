process ARRIBA_DOWNLOAD {
    tag "arriba"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gnu-wget=1.18" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h5bf99c6_5"
    } else {
        container "quay.io/biocontainers/gnu-wget:1.18--h5bf99c6_5"
    }

    output:
    path "versions.yml"   , emit: versions
    path "*"              , emit: reference

    script:
    """
    wget https://github.com/suhrig/arriba/releases/download/v2.3.0/arriba_v2.3.0.tar.gz -O arriba_v2.3.0.tar.gz
    tar -xzvf arriba_v2.3.0.tar.gz
    rm arriba_v2.3.0.tar.gz
    mv arriba_v2.3.0/database/* .
    rm -r arriba_v2.3.0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3 > versions.yml)
    END_VERSIONS
    """

    stub:
    """
    touch blacklist_hg38_GRCh38_v2.3.0.tsv.gz
    touch protein_domains_hg38_GRCh38_v2.3.0.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3 > versions.yml)
    END_VERSIONS
    """
}
