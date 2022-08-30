
process SQUID_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::squid=1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/nfcore/rnafusion:squid_1.5-star2.7.1a' :
        'docker.io/nfcore/rnafusion:squid_1.5-star2.7.1a' }"



    input:
    tuple val(meta), path(txt)
    path gtf

    output:
    tuple val(meta), path("*annotated.txt")    , emit: fusions_annotated
    path  "versions.yml"                        , emit: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    AnnotateSQUIDOutput.py $gtf $txt ${prefix}.squid.fusions.annotated.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        squid: \$(echo \$(squid --version 2>&1) | sed 's/v//')
    END_VERSIONS
    """
}
