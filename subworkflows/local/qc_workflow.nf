//
// Check input samplesheet and get read channels
//

include { MOSDEPTH }                                   from '../../modules/nf-core/mosdepth/main'
include { PICARD_COLLECTRNASEQMETRICS }                from '../../modules/local/picard/collectrnaseqmetrics/main'
include { PICARD_MARKDUPLICATES }                      from '../../modules/nf-core/picard/markduplicates/main'

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

        ch_bam_sorted_indexed
            .map { meta, bam, bai ->
            return [meta, bam, bai, []]
            }.set { ch_bam_sorted_indexed_bed }

        MOSDEPTH(ch_bam_sorted_indexed_bed, ch_fasta)
        ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
        ch_mosdepth_summary = Channel.empty().mix(MOSDEPTH.out.summary_txt)
        ch_mosdepth_global = Channel.empty().mix(MOSDEPTH.out.global_txt)

        PICARD_COLLECTRNASEQMETRICS(ch_bam_sorted_indexed, ch_refflat, ch_rrna_interval)
        ch_versions = ch_versions.mix(PICARD_COLLECTRNASEQMETRICS.out.versions)
        ch_rnaseq_metrics = Channel.empty().mix(PICARD_COLLECTRNASEQMETRICS.out.metrics)

        PICARD_MARKDUPLICATES(ch_bam_sorted, ch_fasta, ch_fai)
        ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)
        ch_duplicate_metrics = Channel.empty().mix(PICARD_MARKDUPLICATES.out.metrics)


    emit:
        versions            = ch_versions.ifEmpty(null)
        mosdepth_summary    = ch_mosdepth_summary
        mosdepth_gobal      = ch_mosdepth_global
        rnaseq_metrics      = ch_rnaseq_metrics
        duplicate_metrics   = ch_duplicate_metrics

}

