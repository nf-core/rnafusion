//
// Check input samplesheet and get read channels
//

include { PICARD_COLLECTRNASEQMETRICS }                from '../../modules/local/picard/collectrnaseqmetrics/main'
include { GATK4_MARKDUPLICATES }                       from '../../modules/nf-core/gatk4/markduplicates/main'
include { PICARD_COLLECTINSERTSIZEMETRICS }            from '../../modules/nf-core/picard/collectinsertsizemetrics/main'

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

        PICARD_COLLECTRNASEQMETRICS(ch_bam_sorted_indexed, ch_refflat, ch_rrna_interval)
        ch_versions = ch_versions.mix(PICARD_COLLECTRNASEQMETRICS.out.versions)
        ch_rnaseq_metrics = Channel.empty().mix(PICARD_COLLECTRNASEQMETRICS.out.metrics)

        GATK4_MARKDUPLICATES(ch_bam_sorted, ch_fasta.map { meta, fasta -> [ fasta ]}, ch_fai.map { meta, fasta_fai -> [ fasta_fai ]})
        ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions)
        ch_duplicate_metrics = Channel.empty().mix(GATK4_MARKDUPLICATES.out.metrics)

        PICARD_COLLECTINSERTSIZEMETRICS(ch_bam_sorted)
        ch_versions = ch_versions.mix(PICARD_COLLECTINSERTSIZEMETRICS.out.versions)
        ch_insertsize_metrics = Channel.empty().mix(PICARD_COLLECTINSERTSIZEMETRICS.out.metrics)


    emit:
        versions            = ch_versions
        rnaseq_metrics      = ch_rnaseq_metrics
        duplicate_metrics   = ch_duplicate_metrics
        insertsize_metrics  = ch_insertsize_metrics

}

