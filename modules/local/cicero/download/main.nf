process CICERO_DOWNLOAD {
    tag "cicero"
    label 'process_low'

    conda "bioconda::gnu-wget=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h5bf99c6_5' :
        'quay.io/biocontainers/gnu-wget:1.18--h5bf99c6_5' }"

    output:
    path "versions.yml"   , emit: versions
    path "*"              , emit: reference

    script:
    """
    wget https://zenodo.org/records/5088371/files/reference.tar.gz?download=1 -O GRCh38_no_alt.tar.gz
    tar -xzvf GRCh38_no_alt.tar.gz
    rm GRCh38_no_alt.tar.gz
    mv reference/Homo_sapiens/GRCh38_no_alt .
    rm -r reference

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3 > versions.yml)
    END_VERSIONS
    """

    stub:
    """
    mkdir GRCh38_no_alt
    touch fake_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3 > versions.yml)
    END_VERSIONS
    """
}
