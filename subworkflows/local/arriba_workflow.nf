//
// Check input samplesheet and get read channels
//

include { STAR_ALIGN as STAR_FOR_ARRIBA }    from '../../modules/nf-core/modules/star/align/main'
include { ARRIBA }                           from '../../modules/nf-core/modules/arriba/main'


workflow ARRIBA_WORKFLOW {
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

        STAR_FOR_ARRIBA( reads, index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center )
        ch_versions = ch_versions.mix(STAR_FOR_ARRIBA.out.versions)

        // reads_align = reads.join(STAR_FOR_STARFUSION.out.junction )

        ARRIBA ( STAR_FOR_ARRIBA.out.bam, fasta, gtf )
        ch_versions = ch_versions.mix(ARRIBA.out.versions)


    emit:
        ARRIBA.out.fusions
        ARRIBA.out.fusions_fail
        versions = ch_versions.ifEmpty(null)
    }

