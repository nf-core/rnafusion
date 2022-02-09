//
// Check input samplesheet and get read channels
//

include { KALLISTO_QUANT    }     from '../../modules/local/kallisto/quant/main'
include { PIZZLY            }     from '../../modules/local/pizzly/detect/main'


workflow PIZZLY_WORKFLOW {
    take:
        reads
        index
        gtf
        transcript

    main:
        ch_versions = Channel.empty()

        KALLISTO_QUANT( reads, index )
        ch_versions = ch_versions.mix(KALLISTO_QUANT.out.versions)



        PIZZLY( KALLISTO_QUANT.out.txt, transcript, gtf )
        // ch_versions = ch_versions.mix(ARRIBA.out.versions)


    emit:
        KALLISTO_QUANT.out.txt
        // ARRIBA.out.fusions_fail
        // versions = ch_versions.ifEmpty(null)
    }

