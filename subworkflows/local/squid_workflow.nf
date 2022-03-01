//
// Check input samplesheet and get read channels
//

include { STAR_ALIGN as STAR_FOR_SQUID }    from '../../modules/nf-core/modules/star/align/main'
// include { ARRIBA }                           from '../../modules/nf-core/modules/arriba/main'


workflow SQUID_WORKFLOW {
    take:
        fasta
        index
        gtf

    main:
        ch_versions = Channel.empty()

        star_ignore_sjdbgtf = false
        seq_platform = false
        seq_center = false

        STAR_FOR_SQUID( fasta, index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center )
        ch_versions = ch_versions.mix(STAR_FOR_ARRIBA.out.versions)

        ARRIBA ( STAR_FOR_ARRIBA.out.bam, fasta, gtf )
        ch_versions = ch_versions.mix(ARRIBA.out.versions)


    emit:
        ARRIBA.out.fusions
        ARRIBA.out.fusions_fail
        versions = ch_versions.ifEmpty(null)
    }

