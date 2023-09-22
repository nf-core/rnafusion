//
// Check input samplesheet and get read channels
//

include { QUALIMAP_RNASEQ }                            from '../../modules/nf-core/qualimap/rnaseq/main'
include { PICARD_COLLECTRNASEQMETRICS }                from '../../modules/local/picard/collectrnaseqmetrics/main'
include { GATK4_MARKDUPLICATES }                      from '../../modules/nf-core/gatk4/markduplicates/main'

workflow QC_WORKFLOW {
    take:
        ch_bam_sorted
        ch_bam_sorted_indexed
        ch_chrgtf
        ch_refflat
        ch_fasta
        ch_fai
        ch_rrna_interval

    main:
        ch_versions = Channel.empty()

        QUALIMAP_RNASEQ(ch_bam_sorted, ch_chrgtf)
        ch_versions = ch_versions.mix(QUALIMAP_RNASEQ.out.versions)
        ch_qualimap_qc = Channel.empty().mix(QUALIMAP_RNASEQ.out.results)

        PICARD_COLLECTRNASEQMETRICS(ch_bam_sorted_indexed, ch_refflat, ch_rrna_interval)
        ch_versions = ch_versions.mix(PICARD_COLLECTRNASEQMETRICS.out.versions)
        ch_rnaseq_metrics = Channel.empty().mix(PICARD_COLLECTRNASEQMETRICS.out.metrics)

        GATK4_MARKDUPLICATES(ch_bam_sorted, ch_fasta, ch_fai)
        ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions)
        ch_duplicate_metrics = Channel.empty().mix(GATK4_MARKDUPLICATES.out.metrics)


    emit:
        versions            = ch_versions.ifEmpty(null)
        qualimap_qc         = ch_qualimap_qc
        rnaseq_metrics      = ch_rnaseq_metrics
        duplicate_metrics   = ch_duplicate_metrics

}

