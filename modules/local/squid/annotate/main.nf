
process SQUID_ANNOTATE {
    tag "squid"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::squid=1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/nfcore/rnafusion:squid_1.5-star2.7.1a' :
        'docker.io/nfcore/rnafusion:squid_1.5-star2.7.1a' }"



    input:
    tuple val(meta), path(txt)
    path gtf

    output:
    tuple val(meta), path("*_annotated.txt"), emit: fusions_annotated
    path  "versions.yml"          , emit: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    AnnotateSQUIDOutput.py $gtf $txt ${prefix}_fusions_annotated.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        squid: \$(echo \$(squid --version 2>&1) | sed 's/v//')
    END_VERSIONS
    """
}
