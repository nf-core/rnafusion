include { AGAT_CONVERTSPGFF2TSV     }                     from '../../modules/nf-core/agat/convertspgff2tsv/main'
include { ARRIBA_VISUALISATION     }                      from '../../modules/local/arriba/visualisation/main'
include { CAT_CAT }                                       from '../../modules/nf-core/cat/cat/main'
include { VCF_COLLECT }                                   from '../../modules/local/vcf_collect/main'
include { FUSIONINSPECTOR     }                           from '../../modules/local/fusioninspector/main'

workflow FUSIONINSPECTOR_WORKFLOW {
    take:
        reads
        fusion_list
        fusion_list_filtered
        fusionreport_out
        fusionreport_csv
        bam_sorted_indexed
        ch_gtf
        ch_arriba_ref_protein_domains
        ch_arriba_ref_cytobands
        ch_hgnc_ref
        ch_hgnc_date

    main:
        ch_versions = Channel.empty()
        ch_arriba_visualisation = Channel.empty()
        index ="${params.starfusion_ref}"

        ch_fusion_list = ( params.tools_cutoff > 1 ? fusion_list_filtered : fusion_list )
        .branch{
            no_fusions: it[1].size() == 0
            fusions: it[1].size() > 0
        }

        if (params.whitelist)  {
            ch_whitelist = ch_fusion_list.fusions.combine(Channel.value(file(params.whitelist, checkIfExists:true)))
                            .map { meta, fusions, whitelist -> [ meta, [fusions, whitelist] ] }

            CAT_CAT(ch_whitelist) // fusioninspector takes care of possible duplicates
            ch_versions = ch_versions.mix(CAT_CAT.out.versions)
            ch_reads_fusion = reads.join(CAT_CAT.out.file_out )
        }
        else {
            ch_reads_fusion = reads.join(ch_fusion_list.fusions )
        }

        FUSIONINSPECTOR( ch_reads_fusion, index)
        ch_versions = ch_versions.mix(FUSIONINSPECTOR.out.versions)

        AGAT_CONVERTSPGFF2TSV(FUSIONINSPECTOR.out.out_gtf)
        ch_versions = ch_versions.mix(AGAT_CONVERTSPGFF2TSV.out.versions)

        fusion_data = FUSIONINSPECTOR.out.tsv_coding_effect.join(AGAT_CONVERTSPGFF2TSV.out.tsv).join(fusionreport_out).join(fusionreport_csv)
        VCF_COLLECT(fusion_data, ch_hgnc_ref, ch_hgnc_date)
        ch_versions = ch_versions.mix(VCF_COLLECT.out.versions)

        if ((params.starfusion || params.all || params.stringtie) && !params.fusioninspector_only && !params.skip_vis) {
            ch_bam_sorted_indexed_fusions = bam_sorted_indexed.join(FUSIONINSPECTOR.out.tsv)
            ARRIBA_VISUALISATION(ch_bam_sorted_indexed_fusions, ch_gtf, ch_arriba_ref_protein_domains, ch_arriba_ref_cytobands)
            ch_versions = ch_versions.mix(ARRIBA_VISUALISATION.out.versions)
            ch_arriba_visualisation = ARRIBA_VISUALISATION.out.pdf
        }

    emit:
        ch_arriba_visualisation
        versions             = ch_versions
}

