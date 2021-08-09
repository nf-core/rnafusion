//
// Uncompress and prepare reference genome files
//

params.genome_options       = [:]
params.star_index_options   = [:]

include { GUNZIP as GUNZIP_FASTA        } from '../../modules/nf-core/modules/gunzip/main'                  addParams( options: params.genome_options )
include { UNTAR as UNTAR_STAR_INDEX     } from '../../modules/nf-core/modules/untar/main'                   addParams( options: params.star_index_options   )
include { SAMTOOLS_FAIDX                } from '../../modules/nf-core/modules/samtools/faidx/main'          addParams( options: params.genome_options )
include { STAR_GENOMEGENERATE           } from '../../modules/nf-core/modules/star/genomegenerate/main'     addParams( options: params.star_index_options )

workflow PREPARE_GENOME {

    take:
    prepare_tool_indices

    main:

    // Uncompress genome fasta file if required
    if (params.fasta.endsWith('.gz'))   ch_fasta = GUNZIP_FASTA ( params.fasta ).gunzipe
    else                                ch_fasta = file(params.fasta)

    // Index the genome fasta
    ch_fasta_fai = Channel.empty()
    if (params.fasta_fai)               ch_fasta_fai = file(params.fasta_fai)
    else                                ch_fasta_fai = SAMTOOLS_FAIDX(ch_fasta).fai

    // Get the GTF file
    if (params.gtf) {
        if (params.gtf.endsWith('.gz')) {
            ch_gtf = GUNZIP_GTF ( params.gtf ).gunzip
        } else {
            ch_gtf = file(params.gtf)
        }
    }

    // Uncompress STAR index or generate from scratch if required
    ch_star_index   = Channel.empty()
    ch_star_version = Channel.empty()
    if ('star' in prepare_tool_indices) {
        if (params.star_index) {
            if (params.star_index.endsWith('.tar.gz')) {
                ch_star_index = UNTAR_STAR_INDEX ( params.star_index ).untar
            } else {
                ch_star_index = file(params.star_index)
            }
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
            if(file("${ch_sf_genome}/AnnotFilterRule.pm").exists()){
                ch_starfusion_resource = file(ch_sf_genome)
            }
        }
        else{
            ch_starfusion_resource = STARFUSION_DOWNLOADGENOME (params.genome).out.reference
        }
    }

    emit:
    fasta               = ch_fasta            // path: genome.fasta
    fai                 = ch_fasta_fai        // path: genome.fasta.fai
    gtf                 = ch_gtf              // path: genome.gtf
    star_index          = ch_star_index       // path: star/index/
    starfusion_resource = ch_starfusion_resource

}



