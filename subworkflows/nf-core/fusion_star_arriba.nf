/*
 * Fusion calling with STAR and Arriba
 */

params.star_align_options = [:]
params.arriba_options = [:]

include { STAR_ALIGN }  from '../../modules/nf-core/modules/star/align/main.nf'    addParams( options: params.star_align_options)
include { ARRIBA }      from '../../modules/nf-core/modules/arriba/main.nf'        addParams( options: params.arriba_options )

workflow FUSION_STAR_ARRIBA {
    take:
    reads           // channel: [ val(meta), [ reads ] ]
    star_index      // channel: STAR index
    genome_fasta    // channel: genome FASTA file
    genome_gtf      // channel: genome GTF file

    main:

    //RNAseq alignment with STAR
    ch_genome_bam = Channel.empty()
    ch_star_version = Channel.empty()
    STAR_ALIGN (
        reads,
        star_index,
        genome_gtf
    )
    ch_genome_bam   = STAR_ALIGN.out.bam
    ch_star_version = STAR_ALIGN.out.version

    //Fusion calling with Arriba
    ch_arriba_fusions        = Channel.empty()
    ch_arriba_failed_fusions = Channel.empty()
    ch_arriba_version        = Channel.empty()
    ARRIBA (
        ch_genome_bam,
        genome_fasta,
        genome_gtf
    )
    ch_arriba_fusions        = ARRIBA.out.fusions
    ch_arriba_failed_fusions = ARRIBA.out.fusions_fail
    ch_arriba_version        = ARRIBA.out.version

    emit:
    bam          = ch_genome_bam            // channel: [ val(meta), [ bam ] ]
    fusions      = ch_arriba_fusions        // channel: [ val(meta), [ fusions ]
    fusions_fail = ch_arriba_failed_fusions // channel: [ val(meta), [ fusions_fail ]
    ch_star_version
    ch_arriba_version
}
