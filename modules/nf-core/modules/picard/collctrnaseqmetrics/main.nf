process PICARD_COLLECTRNASEQMETRICS{

    input:
        tuple val(sample), path(bam), path(bai)
        path(refflat)
        path(rrna_intervals)

    output:
        tuple val(sample), path("${sample}_rna_metrics.txt"), emit metrics

    script:
    def strandedness = ''
    if (params.strandedness == 'forward') {
        strandedness = '--STRAND_SPECIFICITY FIRST_READ_TRANSCRIPTION_STRAND'
    } else if (params.strandedness == 'reverse') {
        strandedness = '--STRAND_SPECIFICITY SECOND_READ_TRANSCRIPTION_STRAND'
    }
    def rrna = rrna_intervals == [] ? '' : "--RIBOSOMAL_INTERVALS ${rrna_intervals}"
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard CollectWgsMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    picard CollectRnaSeqMetrics \\
        --TMP_DIR ./tmp \\
        ${strandedness} \\
        ${rrna} \\
        --REF_FLAT ${refflat} \\
        --INPUT ${bam} \\
        --OUTPUT ${prefix}_rna_metrics.txt \\
    """
}