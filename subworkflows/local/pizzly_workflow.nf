//
// Check input samplesheet and get read channels
//

include { KALLISTO_QUANT }     from '../../modules/local/kallisto/quant/main'
// include { ARRIBA }                           from '../../modules/nf-core/modules/arriba/main'


workflow PIZZLY_WORKFLOW {
    take:
        fasta
        index
        transcript
        gtf

    main:
        ch_versions = Channel.empty()

        KALLISTO_QUANT( fasta, index, gtf,)
        // ch_versions = ch_versions.mix(STAR_FOR_ARRIBA.out.versions)

        // ARRIBA ( STAR_F/OR_ARRIBA.out.bam, fasta, gtf )
        // ch_versions = ch_versions.mix(ARRIBA.out.versions)


    emit:
        KALLISTO_QUANT.out.txt
        // ARRIBA.out.fusions_fail
        // versions = ch_versions.ifEmpty(null)
    }

