//
// Check input samplesheet and get read channels
//

include { STAR_ALIGN as STAR_FOR_STARFUSION }    from '../../modules/nf-core/modules/star/align/main'
include { STARFUSION }                           from '../../modules/local/starfusion/detect/main'


workflow STARFUSION_WORKFLOW {
    take:
        reads
        index
        gtf
        starfusion_ref

    main:
        ch_versions = Channel.empty()

        star_ignore_sjdbgtf = false
        seq_platform = false
        seq_center = false

        STAR_FOR_STARFUSION( reads, index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center )
        // ch_versions = ch_versions.mix(STAR_FOR_STARFUSION.out.versions)
        reads_junction = reads.join(STAR_FOR_STARFUSION.out.junction )

        STARFUSION( reads_junction, starfusion_ref)
        // ch_versions = ch_versions.mix(STARFUSION.out.versions, starfusion_ref )

    emit:
        versions        = ch_versions.ifEmpty(null) // channel: [ versions.yml ]

    // versions = ch_versions.ifEmpty(null)
    // STARFUSION.out.fusions
    // STARFUSION.out.abridged
    // STARFUSION.out.coding_effect

}

