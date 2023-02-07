include { FUSIONINSPECTOR     }                           from '../../modules/local/fusioninspector/main'

workflow FUSIONINSPECTOR_WORKFLOW {
    take:
        reads
        fusion_list
        fusion_list_filtered

    main:
        ch_versions = Channel.empty()
        index ="${params.starfusion_ref}"
        ch_fusion_list = params.fusioninspector_filter ? fusion_list_filtered : fusion_list
        reads_fusion = reads.join(ch_fusion_list )

        FUSIONINSPECTOR( reads_fusion, index)
        ch_versions = ch_versions.mix(FUSIONINSPECTOR.out.versions)

    emit:
        versions        = ch_versions.ifEmpty(null)
}

