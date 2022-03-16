//
// Check input samplesheet and get read channels
//

include { STAR_ALIGN as STAR_FOR_STARFUSION }    from '../../modules/nf-core/modules/star/align/main'
include { STARFUSION }                           from '../../modules/local/starfusion/detect/main'
include { GET_PATH }                         from '../../modules/local/getpath/main'


workflow STARFUSION_WORKFLOW {
    take:
        reads

    main:
        ch_versions = Channel.empty()
        ch_dummy_file = file("$baseDir/assets/dummy_file_starfusion.txt", checkIfExists: true)

        if (params.starfusion){
            if (params.starfusion_fusions){
                ch_starfusion_fusions = params.starfusion_fusions
            } else {
                gtf ="${params.ensembl_ref}/Homo_sapiens.GRCh38.${params.ensembl_version}.chr.gtf"
                ref ="${params.starfusion_ref}/ctat_genome_lib_build_dir"

                star_ignore_sjdbgtf = false
                seq_platform = false
                seq_center = false

                STAR_FOR_STARFUSION( reads, params.starindex_ref, gtf, star_ignore_sjdbgtf, seq_platform, seq_center )
                ch_versions = ch_versions.mix(STAR_FOR_STARFUSION.out.versions)
                reads_junction = reads.join(STAR_FOR_STARFUSION.out.junction )

                STARFUSION( reads_junction, ref)
                ch_versions = ch_versions.mix(STARFUSION.out.versions)

                GET_PATH(STARFUSION.out.fusions)
                ch_starfusion_fusions = GET_PATH.out.file
            }
        }
        else {
            ch_starfusion_fusions = ch_dummy_file
        }
    emit:
        fusions         = ch_starfusion_fusions
        versions        = ch_versions.ifEmpty(null)


}

