include { FUSIONINSPECTOR }                           from '../../modules/local/fusioninspector/main'



workflow FUSIONINSPECTOR_ONLY_WORKFLOW {
    take:
        reads


    main:
        ch_versions = Channel.empty()
        index ="${params.starfusion_ref}"
        ch_fusion_list = "${params.fusioninspector_list}"

        FUSIONINSPECTOR( reads, ch_fusion_list , index )
        ch_versions = ch_versions.mix(FUSIONINSPECTOR.out.versions)

    emit:
        versions        = ch_versions.ifEmpty(null)
}

