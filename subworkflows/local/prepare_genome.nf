//
// Uncompress and prepare reference genome files
//

params.genome_options                   = [:]
params.star_index_options               = [:]
params.starfusion_untar_options         = [:]
params.starfusion_download_options      = [:]
params.fusioncatcher_download_options   = [:]

include { GUNZIP as GUNZIP_FASTA                } from '../../modules/nf-core/modules/gunzip/main'                  addParams( options: params.genome_options )
include { UNTAR as UNTAR_STAR_INDEX             } from '../../modules/nf-core/modules/untar/main'                   addParams( options: params.star_index_options )
include { UNTAR as UNTAR_STARFUSION_GENOME      } from '../../modules/nf-core/modules/untar/main'                   addParams( options: params.starfusion_untar_options )
include { SAMTOOLS_FAIDX                        } from '../../modules/nf-core/modules/samtools/faidx/main'          addParams( options: params.genome_options )
include { STAR_GENOMEGENERATE                   } from '../../modules/nf-core/modules/star/genomegenerate/main'     addParams( options: params.star_index_options )
include { STARFUSION_DOWNLOADGENOME             } from '../../modules/local/starfusion/download/main'               addParams( options: params.starfusion_download_options)
include { FUSIONCATCHER_DOWNLOADGENOME          } from '../../modules/local/fusioncatcher/download/main'            addParams( options: params.fusioncatcher_download_options)

workflow PREPARE_GENOME {

    take:
    prepare_tool_indices

    main:

    // Uncompress genome fasta file if required
    ch_fasta = params.fasta.endsWith('.gz') ? GUNZIP_FASTA ( params.fasta ).gunzipe : file(params.fasta)

    // Index the genome fasta
    ch_fasta_fai = Channel.empty()
    ch_fasta_fai = params.fasta_fai ? file(params.fasta_fai) : SAMTOOLS_FAIDX(ch_fasta).fai

    // Get the GTF file
    if (params.gtf) {
        ch_gtf = params.gtf.endsWith('.gz') ? GUNZIP_GTF ( params.gtf ).gunzip : file(params.gtf)
    }

    // Uncompress STAR index or generate from scratch if required
    ch_star_index   = Channel.empty()
    ch_star_version = Channel.empty()
    if ('star' in prepare_tool_indices) {
        if (params.star_index) {
            ch_star_index = params.star_index.endsWith('.tar.gz') ? UNTAR_STAR_INDEX ( params.star_index ).untar : file(params.star_index)
        } else {
            ch_star_index   = STAR_GENOMEGENERATE ( ch_fasta, ch_gtf ).index
            ch_star_version = STAR_GENOMEGENERATE.out.version
        }
    }

    //Check genome data required for STAR-Fusion
    ch_starfusion_resource = Channel.empty()
    if (params.starfusion){
        if (params.starfusion_genome) {
            ch_sf_genome = params.starfusion_genome
            if (ch_sf_genome.endsWith('.tar.gz')){
                root_path = UNTAR_STARFUSION_GENOME (ch_sf_genome).untar
                ch_starfusion_resource = file("${root_path}/ctat_genome_lib_build_dir")
            }
            else {
                ch_starfusion_resource = file(ch_sf_genome)
            }
        }
        else{
            ch_starfusion_resource = STARFUSION_DOWNLOADGENOME (params.genome).out.reference
        }
    }

    //Check genome data required for FusionCatcher
    ch_fusioncatcher_resource = Channel.empty()
    if (params.fusioncatcher){
        if (params.fusioncatcher_genome){
            ch_fusioncatcher_resource = file(params.fusioncatcher_genome)
        }
        else{
            ch_fusioncatcher_resource = FUSIONCATCHER_DOWNLOADGENOME (params.organism).out.reference
        }
    }

    emit:
    fasta                   = ch_fasta                  // path: genome.fasta
    fai                     = ch_fasta_fai              // path: genome.fasta.fai
    gtf                     = ch_gtf                    // path: genome.gtf
    star_index              = ch_star_index             // path: star/index/
    starfusion_resource     = ch_starfusion_resource    // path: starfusion genome
    fusioncatcher_resource  = ch_fusioncatcher_resource // path: fusioncatcher genome


}



