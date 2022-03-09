//
// Check input samplesheet and get read channels
//

include { GET_PATH }                                    from '../../modules/local/getpath/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_FOR_SQUID }    from '../../modules/nf-core/modules/samtools/sort/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_FOR_SQUID }    from '../../modules/nf-core/modules/samtools/view/main'
include { SQUID }                                       from '../../modules/local/squid/detect/main'
include { SQUID_ANNOTATE }                              from '../../modules/local/squid/annotate/main'
include { STAR_ALIGN as STAR_FOR_SQUID }                from '../../modules/local/star/align/main'

workflow SQUID_WORKFLOW {
    take:
        reads
        fast

    main:
        ch_versions = Channel.empty()
        ch_dummy_file = file("$baseDir/assets/dummy_file_squid.txt", checkIfExists: true)

        if (params.squid) {
            if (params.squid_fusions){
                ch_squid_fusions = params.squid_fusions
            } else {
            gtf ="${params.ensembl_ref}/Homo_sapiens.GRCh38.${params.ensembl_version}.gtf"
            star_ignore_sjdbgtf = false
            seq_platform = false
            seq_center = false
            index = params.starindex_ref

            STAR_FOR_SQUID( reads, index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center )
            ch_versions = ch_versions.mix(STAR_FOR_SQUID.out.versions )

            SAMTOOLS_VIEW_FOR_SQUID ( STAR_FOR_SQUID.out.sam, [] )
            ch_versions = ch_versions.mix(SAMTOOLS_VIEW_FOR_SQUID.out.versions )

            SAMTOOLS_SORT_FOR_SQUID ( SAMTOOLS_VIEW_FOR_SQUID.out.bam )
            ch_versions = ch_versions.mix(SAMTOOLS_SORT_FOR_SQUID.out.versions )

            bam_sorted = STAR_FOR_SQUID.out.bam_sorted.join(SAMTOOLS_SORT_FOR_SQUID.out.bam )

            SQUID ( bam_sorted )
            ch_versions = ch_versions.mix(SQUID.out.versions)

            SQUID_ANNOTATE ( SQUID.out.fusions, gtf )
            ch_versions = ch_versions.mix(SQUID_ANNOTATE.out.versions)

            GET_PATH(SQUID_ANNOTATE.out.fusions_annotated)
            ch_squid_fusions = GET_PATH.out.file
            }
        }
        else {
            ch_squid_fusions = ch_dummy_file
        }

    emit:
        fusions  = ch_squid_fusions
        versions = ch_versions.ifEmpty(null)
    }

