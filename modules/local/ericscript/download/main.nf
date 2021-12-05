// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ERICSCRIPT_DOWNLOAD {
    tag "eriscript"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::awscli" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.7"
    } else {
        container "quay.io/bitnami/python:3.7"
    }

    output:
    path "*.version.txt" , emit: version
    path "homo_sapiens/*", emit: reference

    script:
    def software = getSoftwareName(task.process)
    """
    pip install awscli
    aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/Ensembl/GRCh38/Sequence/ericscript_db_homosapiens_ensembl84.tar.bz2 .
    tar jxf ericscript_db_homosapiens_ensembl84.tar.bz2 --strip-components=2
    rm ericscript_db_homosapiens_ensembl84.tar.bz2

    echo \$(aws --version 2>&1) cut -d "/" -f 2 | cut -d " " -f 1 > ${software}.version.txt
    """
}