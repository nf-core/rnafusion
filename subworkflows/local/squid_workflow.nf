include { GET_META }                                    from '../../modules/local/getmeta/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FOR_SQUID}   from '../../modules/nf-core/modules/samtools/index/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_FOR_SQUID }    from '../../modules/nf-core/modules/samtools/sort/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_FOR_SQUID }    from '../../modules/nf-core/modules/samtools/view/main'
include { SQUID }                                       from '../../modules/local/squid/detect/main'
include { SQUID_ANNOTATE }                              from '../../modules/local/squid/annotate/main'
include { STAR_ALIGN as STAR_FOR_SQUID }                from '../../modules/nf-core/modules/star/align/main'

workflow SQUID_WORKFLOW {

    take:
        reads
        ch_gtf
        ch_starindex_ensembl_ref

    main:
        ch_versions = Channel.empty()
        ch_dummy_file = file("$baseDir/assets/dummy_file_squid.txt", checkIfExists: true)

        if (params.squid || params.all) {
            if (params.squid_fusions){
                ch_squid_fusions = GET_META(reads, params.squid_fusions)
            } else {

            STAR_FOR_SQUID( reads, ch_starindex_ensembl_ref, ch_gtf, params.star_ignore_sjdbgtf, params.seq_platform, params.seq_center )
            ch_versions = ch_versions.mix(STAR_FOR_SQUID.out.versions )

            STAR_FOR_SQUID.out.sam
            .map { meta, sam ->
            return [meta, sam, []]
            }.set { sam_indexed }

            SAMTOOLS_VIEW_FOR_SQUID ( sam_indexed, [] )
            ch_versions = ch_versions.mix(SAMTOOLS_VIEW_FOR_SQUID.out.versions )

            SAMTOOLS_SORT_FOR_SQUID ( SAMTOOLS_VIEW_FOR_SQUID.out.bam )
            ch_versions = ch_versions.mix(SAMTOOLS_SORT_FOR_SQUID.out.versions )

            bam_sorted = STAR_FOR_SQUID.out.bam_sorted.join(SAMTOOLS_SORT_FOR_SQUID.out.bam )

            SQUID ( bam_sorted )
            ch_versions = ch_versions.mix(SQUID.out.versions)

            SQUID_ANNOTATE ( SQUID.out.fusions, ch_gtf )
            ch_versions = ch_versions.mix(SQUID_ANNOTATE.out.versions)

            ch_squid_fusions = SQUID_ANNOTATE.out.fusions_annotated
            }
        }
        else {
            ch_squid_fusions = GET_META(reads, ch_dummy_file)
        }

    emit:
        fusions  = ch_squid_fusions
        versions = ch_versions.ifEmpty(null)
    }

