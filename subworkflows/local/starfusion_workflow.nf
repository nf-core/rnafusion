include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FOR_STARFUSION }      from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FOR_STARFUSION_CRAM } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_FOR_STARFUSION }        from '../../modules/nf-core/samtools/view/main'
include { STAR_ALIGN as STAR_FOR_STARFUSION }                    from '../../modules/nf-core/star/align/main'
include { STARFUSION }                                           from '../../modules/local/starfusion/detect/main'

workflow STARFUSION_WORKFLOW {
    take:
        reads
        ch_chrgtf
        ch_starindex_ref
        ch_fasta

    main:
        ch_versions = Channel.empty()
        ch_align = Channel.empty()
        bam_sorted_indexed = Channel.empty()

        ch_dummy_file = file("$baseDir/assets/dummy_file_starfusion.txt", checkIfExists: true)

        if ((params.starfusion || params.all || params.stringtie) && !params.fusioninspector_only) {
            if (params.starfusion_fusions){
                ch_starfusion_fusions = reads.combine(Channel.value(file(params.starfusion_fusions, checkIfExists:true)))
                                        .map { meta, reads, fusions -> [ meta, fusions ] }
            } else {
                STAR_FOR_STARFUSION( reads, ch_starindex_ref, ch_chrgtf, params.star_ignore_sjdbgtf, '', params.seq_center ?: '')
                ch_versions = ch_versions.mix(STAR_FOR_STARFUSION.out.versions)
                ch_align = STAR_FOR_STARFUSION.out.bam_sorted

                SAMTOOLS_INDEX_FOR_STARFUSION(STAR_FOR_STARFUSION.out.bam_sorted)
                ch_versions = ch_versions.mix(SAMTOOLS_INDEX_FOR_STARFUSION.out.versions)
                bam_sorted_indexed = STAR_FOR_STARFUSION.out.bam_sorted.join(SAMTOOLS_INDEX_FOR_STARFUSION.out.bai)
                reads_junction = reads.join(STAR_FOR_STARFUSION.out.junction )

                if (params.cram.contains('starfusion')){
                    SAMTOOLS_VIEW_FOR_STARFUSION (bam_sorted_indexed, ch_fasta, [] )
                    ch_versions = ch_versions.mix(SAMTOOLS_VIEW_FOR_STARFUSION.out.versions)

                    SAMTOOLS_INDEX_FOR_STARFUSION_CRAM (SAMTOOLS_VIEW_FOR_STARFUSION.out.cram)
                    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_FOR_STARFUSION_CRAM.out.versions)
                }
                if (params.starfusion || params.all){
                    STARFUSION( reads_junction, params.starfusion_ref)
                    ch_versions = ch_versions.mix(STARFUSION.out.versions)
                    ch_starfusion_fusions = STARFUSION.out.fusions
                }

                ch_star_stats = STAR_FOR_STARFUSION.out.log_final
                ch_star_gene_count = STAR_FOR_STARFUSION.out.read_per_gene_tab
            }
        }
        else {
            ch_starfusion_fusions = reads.combine(Channel.value(file(ch_dummy_file, checkIfExists:true)))
                                    .map { meta, reads, fusions -> [ meta, fusions ] }
            ch_star_stats = Channel.empty()
            ch_star_gene_count = Channel.empty()
        }
    emit:
        fusions               = ch_starfusion_fusions
        star_stats            = ch_star_stats
        star_gene_count       = ch_star_gene_count
        ch_bam_sorted         = ch_align.ifEmpty([[],[]])
        ch_bam_sorted_indexed = bam_sorted_indexed.ifEmpty([[],[],[]])
        versions              = ch_versions

    }

