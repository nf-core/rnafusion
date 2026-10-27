process HGNC_DOWNLOAD {
    tag "hgnc"
    label 'process_low'

    conda "bioconda::gnu-wget=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3b/3b54fa9135194c72a18d00db6b399c03248103f87e43ca75e4b50d61179994b3/data' :
        'community.wave.seqera.io/library/wget:1.21.4--8b0fcde81c17be5e' }"

    output:
    path "hgnc_complete_set.txt"        , emit: hgnc_ref
    path "HGNC-DB-timestamp.txt"        , emit: hgnc_date
    tuple val("${task.process}"), val('wget'), eval("wget --version | head -1 | cut -d ' ' -f 3"), topic: versions, emit: versions_wget


    script:
    """
    wget --no-check-certificate https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
    date +%Y-%m-%d/%H:%M  > HGNC-DB-timestamp.txt
    """

    stub:
    """
    touch "hgnc_complete_set.txt"
    touch "HGNC-DB-timestamp.txt"
    """

}
