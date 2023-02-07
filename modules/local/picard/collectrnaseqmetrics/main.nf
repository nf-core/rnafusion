process PICARD_COLLECTRNASEQMETRICS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::picard=2.26.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.10--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.10--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(refflat)
    path(rrna_intervals)

    output:
    tuple val(meta), path("*rna_metrics.txt")    , emit: metrics
    path  "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def strandedness = ''
    // def strandedness = '--STRAND_SPECIFICITY FIRST_READ_TRANSCRIPTION_STRAND'
    if ("${meta.strandedness}" == 'forward') {
        strandedness = '--STRAND_SPECIFICITY FIRST_READ_TRANSCRIPTION_STRAND'
    } else if ("${meta.strandedness}" == 'reverse') {
        strandedness = '--STRAND_SPECIFICITY SECOND_READ_TRANSCRIPTION_STRAND'
    } else {
        strandedness = '--STRAND_SPECIFICITY NONE'
    }

    def rrna = rrna_intervals == [] ? '' : "--RIBOSOMAL_INTERVALS ${rrna_intervals}"
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard CollectRnaMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    picard \\
        -Xmx${avail_mem}g \\
        CollectRnaSeqMetrics \\
        --TMP_DIR ./tmp \\
        ${strandedness} \\
        ${rrna} \\
        --REF_FLAT ${refflat} \\
        --INPUT ${bam} \\
        --OUTPUT ${prefix}_rna_metrics.txt \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectRnaMetrics --version 2>&1 | grep -o 'Version.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
