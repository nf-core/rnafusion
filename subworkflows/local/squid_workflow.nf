//
// Check input samplesheet and get read channels
//

include { STAR_ALIGN as STAR_FOR_SQUID }    from '../../modules/local/star/align/main'
include { SAMTOOLS_VIEW }                   from '../../modules/nf-core/modules/samtools/view/main'
include { SAMTOOLS_SORT }                   from '../../modules/nf-core/modules/samtools/sort/main'
include { SQUID }                           from '../../modules/local/squid/detect/main'
include { SQUID_ANNOTATE }                  from '../../modules/local/squid/annotate/main'


workflow SQUID_WORKFLOW {
    take:
        reads
        fasta
        index
        gtf

    main:
        ch_versions = Channel.empty()

        star_ignore_sjdbgtf = false
        seq_platform = false
        seq_center = false

        STAR_FOR_SQUID( reads, index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center )
        ch_versions = ch_versions.mix(STAR_FOR_SQUID.out.versions )

        SAMTOOLS_VIEW ( STAR_FOR_SQUID.out.sam, [] )
        ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions )

        SAMTOOLS_SORT ( SAMTOOLS_VIEW.out.bam )
        ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions )

        bam_sorted = STAR_FOR_SQUID.out.bam_sorted.join(SAMTOOLS_SORT.out.bam )

        SQUID ( bam_sorted )
        ch_versions = ch_versions.mix(SQUID.out.versions)

        SQUID_ANNOTATE ( SQUID.out.fusions, gtf )
        ch_versions = ch_versions.mix(SQUID_ANNOTATE.out.versions)


    emit:
        // ARRIBA.out.fusions
        // ARRIBA.out.fusions_fail
        versions = ch_versions.ifEmpty(null)
    }

