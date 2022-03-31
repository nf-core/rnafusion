//
// Check input samplesheet and get read channels
//

include { FUSIONINSPECTOR }                           from '../../modules/local/fusioninspector/main'


workflow FUSIONINSPECTOR_WORKFLOW {
    take:
        reads
        fusion_list

    main:
        ch_versions = Channel.empty()
        index ="${params.starfusion_ref}"

        FUSIONINSPECTOR( reads, fusion_list , index )
        ch_versions = ch_versions.mix(FUSIONINSPECTOR.out.versions)

    emit:
        versions        = ch_versions.ifEmpty(null)
}

