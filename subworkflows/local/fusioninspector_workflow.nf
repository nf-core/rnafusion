include { ARRIBA_VISUALISATION     }                      from '../../modules/local/arriba/visualisation/main'
include { CAT_CAT }                                       from '../../modules/nf-core/cat/cat/main'
include { MEGAFUSION }                                    from '../../modules/local/megafusion/main'
include { FUSIONINSPECTOR     }                           from '../../modules/local/fusioninspector/main'

workflow FUSIONINSPECTOR_WORKFLOW {
    take:
        reads
        fusion_list
        fusion_list_filtered
        report
        bam_sorted_indexed
        ch_gtf

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

        fusion_data = FUSIONINSPECTOR.out.tsv.join(report)
        MEGAFUSION(fusion_data)

        if ((params.starfusion || params.all || params.stringtie) && !params.fusioninspector_only) {
            bam_sorted_indexed_fusions = bam_sorted_indexed.join(FUSIONINSPECTOR.out.fusions)
            ARRIBA_VISUALISATION(bam_sorted_indexed_fusions, ch_gtf, params.arriba_ref_protein_domain, params.arriba_ref_cytobands)
            ch_versions = ch_versions.mix(ARRIBA_VISUALISATION.out.versions)
        }

    emit:
        versions        = ch_versions.ifEmpty(null)
}

