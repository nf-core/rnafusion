include { FUSIONINSPECTOR     }                           from '../../modules/local/fusioninspector/main'
include { CAT_CAT }                                       from '../../modules/nf-core/cat/cat/main'
include { MEGAFUSION }                                    from '../../modules/local/megafusion/main'

workflow FUSIONINSPECTOR_WORKFLOW {
    take:
        reads
        fusion_list
        fusion_list_filtered

    main:
        ch_versions = Channel.empty()
        index ="${params.starfusion_ref}"
        ch_fusion_list = params.fusioninspector_filter ? fusion_list_filtered : fusion_list

        if (params.whitelist)  {
            ch_whitelist = ch_fusion_list.combine(Channel.value(file(params.whitelist, checkIfExists:true)))
                            .map { meta, fusions, whitelist -> [ meta, [fusions, whitelist] ] }

            CAT_CAT(ch_whitelist) // fusioninspector takes care of possible duplicates
            ch_versions = ch_versions.mix(CAT_CAT.out.versions)

            ch_fusion_list = CAT_CAT.out.file_out
        }

        reads_fusion = reads.join(ch_fusion_list )

        FUSIONINSPECTOR( reads_fusion, index)
        ch_versions = ch_versions.mix(FUSIONINSPECTOR.out.versions)

        MEGAFUSION(FUSIONINSPECTOR.out.tsv)


    emit:
        versions        = ch_versions.ifEmpty(null)
}

