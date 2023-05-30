include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FOR_SQUID }             from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FOR_SQUID_CHIMERIC }    from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_FOR_SQUID_CHIMERIC }      from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_FOR_SQUID_CHIMERIC }      from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_FOR_SQUID_CRAM }          from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_FOR_SQUID_CRAM_CHIMERIC } from '../../modules/nf-core/samtools/view/main'
include { SQUID }                                                  from '../../modules/local/squid/detect/main'
include { SQUID_ANNOTATE }                                         from '../../modules/local/squid/annotate/main'
include { STAR_ALIGN as STAR_FOR_SQUID }                           from '../../modules/nf-core/star/align/main'

workflow SQUID_WORKFLOW {

    take:
        reads
        ch_gtf
        ch_starindex_ensembl_ref
        ch_fasta

    main:
        ch_versions = Channel.empty()
        ch_dummy_file = file("$baseDir/assets/dummy_file_squid.txt", checkIfExists: true)

        if ((params.squid || params.all) && !params.fusioninspector_only) {
            if (params.squid_fusions){
                ch_squid_fusions = reads.combine(Channel.value(file(params.squid_fusions, checkIfExists:true)))
                                    .map { meta, reads, fusions -> [ meta, fusions ] }
            } else {

            STAR_FOR_SQUID(reads, ch_starindex_ensembl_ref, ch_gtf, params.star_ignore_sjdbgtf, '', params.seq_center ?: '')
            ch_versions = ch_versions.mix(STAR_FOR_SQUID.out.versions)

            STAR_FOR_SQUID.out.sam
            .map { meta, sam ->
            return [meta, sam, []]
            }.set { chimeric_sam }


            ch_fasta
            .map { it -> tuple(id:it.baseName, it) }
            .set { ch_fasta_w_meta }

            SAMTOOLS_VIEW_FOR_SQUID_CHIMERIC (chimeric_sam, ch_fasta_w_meta, [])
            ch_versions = ch_versions.mix(SAMTOOLS_VIEW_FOR_SQUID_CHIMERIC.out.versions)

            SAMTOOLS_SORT_FOR_SQUID_CHIMERIC (SAMTOOLS_VIEW_FOR_SQUID_CHIMERIC.out.bam)
            ch_versions = ch_versions.mix(SAMTOOLS_SORT_FOR_SQUID_CHIMERIC.out.versions)

            bam_chimeric = STAR_FOR_SQUID.out.bam_sorted.join(SAMTOOLS_SORT_FOR_SQUID_CHIMERIC.out.bam)

            if (params.cram.contains('squid')){
                SAMTOOLS_INDEX_FOR_SQUID(STAR_FOR_SQUID.out.bam_sorted)
                ch_versions = ch_versions.mix(SAMTOOLS_INDEX_FOR_SQUID.out.versions)
                SAMTOOLS_INDEX_FOR_SQUID_CHIMERIC(SAMTOOLS_SORT_FOR_SQUID_CHIMERIC.out.bam)
                ch_versions = ch_versions.mix(SAMTOOLS_INDEX_FOR_SQUID_CHIMERIC.out.versions)

                bam_sorted_indexed = STAR_FOR_SQUID.out.bam_sorted.join(SAMTOOLS_INDEX_FOR_SQUID.out.bai)
                chimeric_sorted_indexed = SAMTOOLS_SORT_FOR_SQUID_CHIMERIC.out.bam.join(SAMTOOLS_INDEX_FOR_SQUID_CHIMERIC.out.bai)

                SAMTOOLS_VIEW_FOR_SQUID_CRAM (bam_sorted_indexed, ch_fasta_w_meta, [])
                ch_versions = ch_versions.mix(SAMTOOLS_VIEW_FOR_SQUID_CRAM.out.versions)
                SAMTOOLS_VIEW_FOR_SQUID_CRAM_CHIMERIC (chimeric_sorted_indexed, ch_fasta_w_meta, [])
                ch_versions = ch_versions.mix(SAMTOOLS_VIEW_FOR_SQUID_CRAM.out.versions)
            }

            SQUID (bam_chimeric)
            ch_versions = ch_versions.mix(SQUID.out.versions)

            SQUID_ANNOTATE (SQUID.out.fusions, ch_gtf)
            ch_versions = ch_versions.mix(SQUID_ANNOTATE.out.versions)

            ch_squid_fusions = SQUID_ANNOTATE.out.fusions_annotated
            }
        }
        else {
            ch_squid_fusions = reads.combine(Channel.value(file(ch_dummy_file, checkIfExists:true)))
                                .map { meta, reads, fusions -> [ meta, fusions ] }
        }

    emit:
        fusions  = ch_squid_fusions
        versions = ch_versions.ifEmpty(null)
    }

