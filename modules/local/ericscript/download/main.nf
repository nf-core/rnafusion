process ERICSCRIPT_DOWNLOAD {
    tag "eriscript"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gnu-wget=1.18" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h5bf99c6_5"
    } else {
        container "quay.io/biocontainers/gnu-wget:1.18--h5bf99c6_5"
    }

    output:
    path "versions.yml"  , emit: versions
    path "homo_sapiens/*", emit: reference

    script:
    """
    wget http://ngi-igenomes.s3.amazonaws.com/igenomes/Homo_sapiens/Ensembl/GRCh38/Sequence/ericscript_db_homosapiens_ensembl84.tar.bz2
    tar jxf ericscript_db_homosapiens_ensembl84.tar.bz2 --strip-components=2
    rm ericscript_db_homosapiens_ensembl84.tar.bz2

    echo \$(wget -V 2>&1) | grep "GNU Wget" | cut -d" " -f3 > versions.yml
    """
}
