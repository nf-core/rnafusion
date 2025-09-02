include { STRINGTIE_STRINGTIE }    from '../../../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_MERGE }        from '../../../modules/nf-core/stringtie/merge/main'


workflow STRINGTIE_WORKFLOW {
    take:
        bam_sorted // channel: [meta, bam]
        ch_chrgtf  // channel: [meta, gtf]

    main:
        ch_versions = Channel.empty()
        ch_stringtie_gtf = Channel.empty()

        STRINGTIE_STRINGTIE(
            bam_sorted,
            ch_chrgtf.map { _meta, gtf -> [ gtf ]}
        )
        ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions)

        STRINGTIE_STRINGTIE
            .out
            .transcript_gtf
            .map { it -> it[1] }
            .set { stringtie_gtf }
        ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions)


        STRINGTIE_MERGE (stringtie_gtf, ch_chrgtf.map { _meta, gtf -> [ gtf ]})
        ch_versions = ch_versions.mix(STRINGTIE_MERGE.out.versions)
        ch_stringtie_gtf = STRINGTIE_MERGE.out.gtf

    emit:
        stringtie_gtf = ch_stringtie_gtf               // channel: [meta, gtf]
        versions      = ch_versions                    // channel: [ path(versions.yml) ]

    }
