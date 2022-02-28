//
// Check input samplesheet and get read channels
//

include { STAR_ALIGN as STAR_FOR_STARFUSION }    from '../../modules/nf-core/modules/star/align/main'
include { STARFUSION }                           from '../../modules/local/starfusion/detect/main'


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
                index ="${params.starfusion_ref}/ctat_genome_lib_build_dir"

                star_ignore_sjdbgtf = false
                seq_platform = false
                seq_center = false

                STAR_FOR_STARFUSION( reads, index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center )
                ch_versions = ch_versions.mix(STAR_FOR_STARFUSION.out.versions)
                reads_junction = reads.join(STAR_FOR_STARFUSION.out.junction )

                STARFUSION( reads_junction, starfusion_ref)
                ch_versions = ch_versions.mix(STARFUSION.out.versions, params.starfusion_ref )

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

