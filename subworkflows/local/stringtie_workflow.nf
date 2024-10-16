include { STRINGTIE_STRINGTIE }    from '../../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_MERGE }        from '../../modules/nf-core/stringtie/merge/main'


workflow STRINGTIE_WORKFLOW {
    take:
        bam_sorted
        ch_gtf

    main:
        ch_versions = Channel.empty()
        ch_stringtie_gtf = Channel.empty()

        if ((params.stringtie || params.all) && !params.fusioninspector_only) {
            STRINGTIE_STRINGTIE(bam_sorted, ch_gtf.map { meta, gtf -> [ gtf ]})
            ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions)

            STRINGTIE_STRINGTIE
                .out
                .transcript_gtf
                .map { it -> it[1] }
                .set { stringtie_gtf }
            ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions)


            STRINGTIE_MERGE (stringtie_gtf, ch_gtf.map { meta, gtf -> [ gtf ]})
            ch_versions = ch_versions.mix(STRINGTIE_MERGE.out.versions)
            ch_stringtie_gtf = STRINGTIE_MERGE.out.gtf
        }

    emit:
        stringtie_gtf = ch_stringtie_gtf.ifEmpty(null)
        versions      = ch_versions

    }

