include { TRIMGALORE                  } from '../../modules/nf-core/modules/trimgalore/main'

workflow TRIM_WORKFLOW {
    take:
        reads

    main:
        ch_versions = Channel.empty()

        if (params.trim) {

            TRIMGALORE (
                reads
            )
            ch_versions = ch_versions.mix(TRIMGALORE.out.versions)
            ch_reads = TRIMGALORE.out.reads
        }
        else {
            ch_reads = reads
        }

    emit:
        ch_trim_reads           = ch_reads
        versions        = ch_versions.ifEmpty(null)
    }

