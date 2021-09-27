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

    conda (params.enable_conda ? "bioconda::ericscript=0.5.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ericscript:0.5.5--hdfd78af_5"
    } else {
        container "quay.io/biocontainers/ericscript:0.5.5--hdfd78af_5"
    }

    output:
    path "*.version.txt" , emit: version
    path "homo_sapiens/*", emit: reference

    script:
    def software = getSoftwareName(task.process)
    """
    wget https://raw.githubusercontent.com/circulosmeos/gdown.pl/dfd6dc910a38a42d550397bb5c2335be2c4bcf54/gdown.pl && chmod +x gdown.pl
    ./gdown.pl "https://drive.google.com/uc?export=download&confirm=qgOc&id=0B9s__vuJPvIiUGt1SnFMZFg4TlE" ericscript_db_homosapiens_ensembl84.tar.bz2
    tar jxf ericscript_db_homosapiens_ensembl84.tar.bz2 --strip-components=2
    rm ./gdown.pl ericscript_db_homosapiens_ensembl84.tar.bz2

    echo \$(wget -V 2>&1) | grep "GNU Wget" | cut -d" " -f3 > ${software}.version.txt
    """
}