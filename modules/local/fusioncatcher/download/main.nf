// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FUSIONCATCHER_DOWNLOAD {
    tag 'fusioncatcher'
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::fusioncatcher=1.33" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/fusioncatcher:1.33--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/fusioncatcher:1.33--hdfd78af_1"
    }

    output:
    path "human_v102"       , emit: reference
    path "*.version.txt"    , emit: version

    script:
    def software = getSoftwareName(task.process)
    def url = "http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.aa"
    """
    if curl --output /dev/null --silent --head --fail $url; then
        wget $url \\
        wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ab
        wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ac
        wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ad
        cat human_v102.tar.gz.* | tar xz
        rm human_v98.tar*
    else
        fusioncatcher-build \\
            -g homo_sapiens \\
            -o human_v102 \\
            $options.args
    fi
    
    fusioncatcher --version | sed 's/fusioncatcher.py //' > ${software}.version.txt
    """
}
