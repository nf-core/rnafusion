include { FUSIONREPORT }                           from '../../modules/local/fusionreport/detect/main'


workflow FUSIONREPORT_WORKFLOW {
    take:
        reads
        fusionreport_ref
        arriba_fusions
        pizzly_fusions
        squid_fusions
        starfusion_fusions
        fusioncatcher_fusions

    main:
        ch_versions = Channel.empty()

        reads_fusions = reads
        .join(arriba_fusions, remainder: true)
        .join(pizzly_fusions, remainder: true)
        .join(squid_fusions, remainder: true)
        .join(starfusion_fusions, remainder: true)
        .join(fusioncatcher_fusions, remainder: true)

        FUSIONREPORT( reads_fusions, fusionreport_ref)
        ch_versions = ch_versions.mix(FUSIONREPORT.out.versions)

    emit:
        versions        = ch_versions.ifEmpty(null)
        fusion_list     = FUSIONREPORT.out.fusion_list

}

