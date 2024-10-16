process HGNC_DOWNLOAD {
    tag "hgnc"
    label 'process_low'

    conda "bioconda::gnu-wget=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h5bf99c6_5' :
        'quay.io/biocontainers/gnu-wget:1.18--h5bf99c6_5' }"

    input:

    output:
    path "hgnc_complete_set.txt"        , emit: hgnc_ref
    path "HGNC-DB-timestamp.txt"        , emit: hgnc_date

    path "versions.yml"   , emit: versions


    script:
    """
    wget https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt
    date +%Y-%m-%d/%H:%M  > HGNC-DB-timestamp.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3 > versions.yml)
    END_VERSIONS
    """

    stub:
    """
    touch "hgnc_complete_set.txt"
    touch "HGNC-DB-timestamp.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3 > versions.yml)
    END_VERSIONS
    """

}
