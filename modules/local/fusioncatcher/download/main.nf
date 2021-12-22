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
    path "human_${human_version}/*"    , emit: reference
    path "*.version.txt"               , emit: version

    script:
    def software = getSoftwareName(task.process)
    def human_version = "v102"
    def url = "http://sourceforge.net/projects/fusioncatcher/files/data/human_${human_version}.tar.gz.aa"
    """
    if wget --spider "$url" 2>/dev/null; then
        wget $options.args $url
        wget $options.args http://sourceforge.net/projects/fusioncatcher/files/data/human_${human_version}.tar.gz.ab
        wget $options.args http://sourceforge.net/projects/fusioncatcher/files/data/human_${human_version}.tar.gz.ac
        wget $options.args http://sourceforge.net/projects/fusioncatcher/files/data/human_${human_version}.tar.gz.ad
        cat human_${human_version}.tar.gz.* | tar xz
        rm human_${human_version}.tar*
    else
        fusioncatcher-build \\
            -g homo_sapiens \\
            -o human_${human_version} \\
            $options.args2
    fi

    fusioncatcher --version | sed 's/fusioncatcher.py //' > ${software}.version.txt
    """
}
