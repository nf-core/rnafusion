process RRNA_TRANSCRIPTS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.12.2"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
        tuple val(meta), path(gtf)

    output:
        tuple val(meta), path("*rrna_intervals.gtf")   , emit: rrna_gtf
        path "versions.yml"                            , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnafusion/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    get_rrna_transcripts.py $gtf ${prefix}_rrna_intervals.gtf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_rrna_intervals.gtf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
