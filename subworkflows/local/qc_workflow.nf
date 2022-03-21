//
// Check input samplesheet and get read channels
//

include { QUALIMAP_RNASEQ }    from '../../modules/nf-core/modules/qualimap/rnaseq/main'
include { GTF_TO_REFFLAT }     from '../../modules/local/uscs/gtftogenepred/main'


workflow QC_WORKFLOW {
    take:
        bam

    main:
        ch_versions = Channel.empty()
        ch_qualimap_qc = Channel.empty()

        if (params.starfusion && params.qc){
            gtf ="${params.ensembl_ref}/Homo_sapiens.GRCh38.${params.ensembl_version}.gtf"

            QUALIMAP_RNASEQ(bam, gtf)
            ch_versions = ch_versions.mix(QUALIMAP_RNASEQ.out.versions)
            ch_qualimap_qc = QUALIMAP_RNASEQ.out.results.ifEmpty(null)

            SAMTOOLS_INDEX_FOR_QC(bam)
            ch_versions = ch_versions.mix(SAMTOOLS_INDEX_FOR_QC.out.versions)

            bam_indexed = bam.join(SAMTOOLS_INDEX_FOR_QC.out.bai)

            GTF_TO_REFFLAT(gtf)
            PICARD_COLLECTRNASEQMETRICS(bam_indexed, GTF_TO_REFFLAT.out.refflat, [])
        }
    emit:
        versions        = ch_versions.ifEmpty(null)
        qualimap_qc     = ch_qualimap_qc.ifEmpty(null)
        rnaseq_metrics  = PICARD_COLLECTRNASEQMETRICS.out.metrics
}

