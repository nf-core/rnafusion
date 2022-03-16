//
// Check input samplesheet and get read channels
//

include { ARRIBA }                                      from '../../modules/nf-core/modules/arriba/main'
include { ARRIBA_VISUALISATION }                        from '../../modules/local/arriba/visualisation/main'
include { GET_PATH as GET_PATH_ARRIBA }                 from '../../modules/local/getpath/main'
include { GET_PATH as GET_PATH_ARRIBA_FAIL }            from '../../modules/local/getpath/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_FOR_ARRIBA }   from '../../modules/nf-core/modules/samtools/sort/main'
include { SAMTOOLS_INDEX }                              from '../../modules/nf-core/modules/samtools/index/main'
include { STAR_ALIGN as STAR_FOR_ARRIBA }               from '../../modules/nf-core/modules/star/align/main'


workflow ARRIBA_WORKFLOW {
    take:
        reads
        fasta
        index

    main:
        ch_versions = Channel.empty()
        ch_dummy_file = file("$baseDir/assets/dummy_file_arriba.txt", checkIfExists: true)

        if (params.arriba) {
            gtf ="${params.ensembl_ref}/Homo_sapiens.GRCh38.${params.ensembl_version}.gtf"
            star_ignore_sjdbgtf = false
            seq_platform = false
            seq_center = false

            STAR_FOR_ARRIBA( reads, index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center )
            ch_versions = ch_versions.mix(STAR_FOR_ARRIBA.out.versions)

            SAMTOOLS_SORT_FOR_ARRIBA(STAR_FOR_ARRIBA.out.bam)
            ch_versions = ch_versions.mix(SAMTOOLS_SORT_FOR_ARRIBA.out.versions)

            SAMTOOLS_INDEX(SAMTOOLS_SORT_FOR_ARRIBA.out.bam)
            ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

            bam_indexed = SAMTOOLS_SORT_FOR_ARRIBA.out.bam.join(SAMTOOLS_INDEX.out.bai)

            if (params.arriba_fusions) {
                ch_arriba_fusions = params.arriba_fusions
                ch_arriba_fusion_fail = ch_dummy_file
            } else {
                ARRIBA ( STAR_FOR_ARRIBA.out.bam, fasta, gtf )
                ch_versions = ch_versions.mix(ARRIBA.out.versions)

                GET_PATH_ARRIBA(ARRIBA.out.fusions)
                ch_arriba_fusions = GET_PATH_ARRIBA.out.file

                GET_PATH_ARRIBA_FAIL(ARRIBA.out.fusions_fail)
                ch_arriba_fusion_fail = GET_PATH_ARRIBA_FAIL.out.file
            }

            ARRIBA_VISUALISATION(bam_indexed, ch_arriba_fusions, params.arriba_ref, gtf)
            ch_versions = ch_versions.mix(ARRIBA_VISUALISATION.out.versions)

            ch_arriba_visualisation = ARRIBA_VISUALISATION.out.pdf

        }
        else {
            ch_arriba_fusions       = ch_dummy_file
            ch_arriba_fusion_fail   = ch_dummy_file
            ch_arriba_visualisation = ch_dummy_file
        }

    emit:
        fusions         = ch_arriba_fusions
        fusions_fail    = ch_arriba_fusion_fail
        versions        = ch_versions.ifEmpty(null)
        pdf             = ch_arriba_visualisation
    }

