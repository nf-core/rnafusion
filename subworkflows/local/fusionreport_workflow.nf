//
// Check input samplesheet and get read channels
//

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

        FUSIONREPORT( reads, fusionreport_ref, arriba_fusions, pizzly_fusions, squid_fusions, starfusion_fusions, fusioncatcher_fusions )
        ch_versions = ch_versions.mix(FUSIONREPORT.out.versions)

    emit:
        versions        = ch_versions.ifEmpty(null)
        fusion_list     = FUSIONREPORT.out.fusion_list

}

