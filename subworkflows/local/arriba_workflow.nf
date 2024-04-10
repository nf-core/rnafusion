include { ARRIBA_ARRIBA }                               from '../../modules/nf-core/arriba/arriba/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FOR_ARRIBA}  from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_FOR_ARRIBA}    from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_FOR_ARRIBA}    from '../../modules/nf-core/samtools/view/main'
include { STAR_ALIGN as STAR_FOR_ARRIBA }               from '../../modules/nf-core/star/align/main'

workflow ARRIBA_WORKFLOW {
    take:
        reads
        ch_gtf
        ch_fasta
        ch_starindex_ref
        ch_arriba_ref_blacklist
        ch_arriba_ref_known_fusions
        ch_arriba_ref_protein_domains

    main:
        ch_versions = Channel.empty()
        ch_dummy_file = file("$baseDir/assets/dummy_file_arriba.txt", checkIfExists: true)

        if ((params.arriba || params.all) && !params.fusioninspector_only) {

            STAR_FOR_ARRIBA( reads, ch_starindex_ref, ch_gtf, params.star_ignore_sjdbgtf, '', params.seq_center ?: '')
            ch_versions = ch_versions.mix(STAR_FOR_ARRIBA.out.versions)

            if (params.arriba_fusions) {
                ch_arriba_fusions = reads.combine( Channel.value( file( params.arriba_fusions, checkIfExists: true ) ) )
                    .map { meta, reads, fusions -> [ meta, fusions ] }
                ch_arriba_fusion_fail = ch_dummy_file
            } else {
                ARRIBA_ARRIBA ( STAR_FOR_ARRIBA.out.bam, ch_fasta, ch_gtf, ch_arriba_ref_blacklist, ch_arriba_ref_known_fusions, [[],[]], [[],[]], ch_arriba_ref_protein_domains )
                ch_versions = ch_versions.mix(ARRIBA_ARRIBA.out.versions)

                ch_arriba_fusions     = ARRIBA_ARRIBA.out.fusions
                ch_arriba_fusion_fail = ARRIBA_ARRIBA.out.fusions_fail.map{ meta, file -> return file}
            }

            if (params.cram.contains('arriba') ){

                SAMTOOLS_SORT_FOR_ARRIBA(STAR_FOR_ARRIBA.out.bam, ch_fasta)
                ch_versions = ch_versions.mix(SAMTOOLS_SORT_FOR_ARRIBA.out.versions )

                SAMTOOLS_VIEW_FOR_ARRIBA(SAMTOOLS_SORT_FOR_ARRIBA.out.bam.map { meta, bam -> [ meta, bam, [] ] }, ch_fasta, [])
                ch_versions = ch_versions.mix(SAMTOOLS_VIEW_FOR_ARRIBA.out.versions )

                SAMTOOLS_INDEX_FOR_ARRIBA(SAMTOOLS_VIEW_FOR_ARRIBA.out.cram)
                ch_versions = ch_versions.mix(SAMTOOLS_INDEX_FOR_ARRIBA.out.versions )

            }

        }
        else {
            ch_arriba_fusions       = reads.combine(Channel.value( file(ch_dummy_file, checkIfExists:true ) ) )
                                        .map { meta, reads, fusions -> [ meta, fusions ] }

            ch_arriba_fusion_fail   = ch_dummy_file
        }

    emit:
        fusions         = ch_arriba_fusions
        fusions_fail    = ch_arriba_fusion_fail
        versions        = ch_versions
    }

