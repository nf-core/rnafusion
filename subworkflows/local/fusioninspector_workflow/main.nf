include { AGAT_CONVERTSPGFF2TSV     }                     from '../../../modules/nf-core/agat/convertspgff2tsv/main'
include { ARRIBA_VISUALISATION     }                      from '../../../modules/nf-core/arriba/visualisation/main'
include { CAT_CAT }                                       from '../../../modules/nf-core/cat/cat/main'
include { VCF_COLLECT }                                   from '../../../modules/local/vcf_collect/main'
include { FUSIONINSPECTOR     }                           from '../../../modules/nf-core/fusioninspector/main'

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
        ch_starfusion_ref
        skip_vis
        skip_vcf
        tools_cutoff
        whitelist

    main:
        ch_versions = Channel.empty()
        ch_arriba_visualisation = Channel.empty()

        ch_fusion_list = ( tools_cutoff > 1 ? fusion_list_filtered : fusion_list )

        if (whitelist)  {
            ch_whitelist = ch_fusion_list.combine(Channel.value(file(whitelist, checkIfExists:true)))
                            .map { meta, fusions, whitelist_file -> [ meta, [fusions, whitelist_file] ] }

            CAT_CAT(ch_whitelist) // fusioninspector takes care of possible duplicates
            ch_versions = ch_versions.mix(CAT_CAT.out.versions)
            ch_reads_fusion = reads.join(CAT_CAT.out.file_out )
        }
        else {
            ch_reads_fusion = reads.join(ch_fusion_list)
        }

        FUSIONINSPECTOR( ch_reads_fusion, ch_starfusion_ref)

        def tsv_nonempty = FUSIONINSPECTOR.out.tsv.filter { _meta, file -> file.exists() && file.size() > 0 }
        def tsv_abridged_nonempty = FUSIONINSPECTOR.out.abridged_tsv.filter { _meta, file -> file.exists() && file.size() > 0 }
        def gtf_nonempty = FUSIONINSPECTOR.out.out_gtf.filter { _meta, file -> file.exists() && file.size() > 0 }
        if (!tsv_nonempty) {
            log.warn("FUSIONINSPECTOR confirmed no fusions, skipping VCF and visualisation steps.")
        }
        if (
            !skip_vcf
        ) {
            AGAT_CONVERTSPGFF2TSV(gtf_nonempty)
            ch_versions = ch_versions.mix(AGAT_CONVERTSPGFF2TSV.out.versions)

            fusion_data = tsv_abridged_nonempty
                .join(AGAT_CONVERTSPGFF2TSV.out.tsv)
                .join(fusionreport_out)
                .join(fusionreport_csv)

            VCF_COLLECT(fusion_data, ch_hgnc_ref, ch_hgnc_date)
        }
        if (
            !skip_vis
        ) {
            ch_bam_sorted_indexed_fusions = bam_sorted_indexed.join(tsv_nonempty)
            ARRIBA_VISUALISATION(
                ch_bam_sorted_indexed_fusions,
                ch_gtf,
                ch_arriba_ref_protein_domains.map { it -> [[id:it.name], it]},
                ch_arriba_ref_cytobands.map { it -> [[id:it.name], it]}
            )
            ch_versions = ch_versions.mix(ARRIBA_VISUALISATION.out.versions)
            ch_arriba_visualisation = ARRIBA_VISUALISATION.out.pdf
        }

    emit:
        ch_arriba_visualisation = ch_arriba_visualisation
        versions             = ch_versions
}
