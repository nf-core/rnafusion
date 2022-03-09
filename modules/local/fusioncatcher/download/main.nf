process FUSIONCATCHER_DOWNLOAD {
    tag "fusioncatcher_download"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::fusioncatcher=1.33" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {

        container "docker.io/clinicalgenomics/fusioncatcher:1.33"
    } else {
        container "docker.io/clinicalgenomics/fusioncatcher:1.33"
    }

    output:
    path "*"                , emit: reference
    path "versions.yml"     , emit: versions

    script:

    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def human_version = "v102"
    def url = "http://sourceforge.net/projects/fusioncatcher/files/data/human_${human_version}.tar.gz.aa"
    """
    if wget --spider "$url" 2>/dev/null; then
        wget $args $url
        wget $args http://sourceforge.net/projects/fusioncatcher/files/data/human_${human_version}.tar.gz.ab
        wget $args http://sourceforge.net/projects/fusioncatcher/files/data/human_${human_version}.tar.gz.ac
        wget $args http://sourceforge.net/projects/fusioncatcher/files/data/human_${human_version}.tar.gz.ad
        cat human_${human_version}.tar.gz.* | tar xz
        rm human_${human_version}.tar*
    else
        fusioncatcher-build \\
            -g homo_sapiens \\
            -o human_${human_version} \\
            $args2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioncatcher: \$(echo \$(fusioncatcher --version 2>&1))
    END_VERSIONS
    """
}
