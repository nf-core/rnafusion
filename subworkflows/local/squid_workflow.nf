//
// Check input samplesheet and get read channels
//

include { STAR_ALIGN as STAR_FOR_SQUID }    from '../../modules/nf-core/modules/star/align/main'
include { SAMTOOLS_SORT }                   from '../../modules/nf-core/modules/samtools/sort/main'
include { SQUID }                           from '../../modules/local/squid/main'


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
        ch_versions = ch_versions.mix(STAR_FOR_SQUID.out.versions)

        SAMTOOLS_SORT ( STAR_FOR_SQUID.out.bam )
        ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

        SQUID ( SAMTOOLS_SORT.out.bam, gtf )
        ch_versions = ch_versions.mix(SQUID.out.versions)


    emit:
        // ARRIBA.out.fusions
        // ARRIBA.out.fusions_fail
        versions = ch_versions.ifEmpty(null)
    }

