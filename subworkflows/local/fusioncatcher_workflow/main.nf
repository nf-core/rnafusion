include { FASTP as FASTP_FOR_FUSIONCATCHER } from '../../../modules/nf-core/fastp/main'
include { FUSIONCATCHER_FUSIONCATCHER }      from '../../../modules/nf-core/fusioncatcher/fusioncatcher/main'

workflow FUSIONCATCHER_WORKFLOW {
    take:
        reads                   // channel [ meta, [ fastqs ] ]
        fusioncatcher_trimming  // boolean
        adapter_fasta           // channel [ path ]
        fusioncatcher_ref       // channel [ meta, path       ]
        fusioncatcher_fusions   // path, string

    main:
        ch_versions   = Channel.empty()

        if (fusioncatcher_fusions){

            ch_fusioncatcher_fusions = reads.combine(Channel.value(file(fusioncatcher_fusions, checkIfExists:true)))
                                        .map { meta, _reads, fusions -> [ meta, fusions ] }
        } else {
            if (fusioncatcher_trimming) {
                FASTP_FOR_FUSIONCATCHER(
                    reads,
                    adapter_fasta,
                    false, // discard_trimmed_pass
                    false, // save_trimmed_fail
                    false  // skip_qc
                )
                ch_versions = ch_versions.mix(FASTP_FOR_FUSIONCATCHER.out.versions)
                reads = FASTP_FOR_FUSIONCATCHER.out.reads
            }

            FUSIONCATCHER_FUSIONCATCHER (
                reads,
                fusioncatcher_ref
            )
            ch_fusioncatcher_fusions = FUSIONCATCHER_FUSIONCATCHER.out.fusions
            ch_versions              = ch_versions.mix(FUSIONCATCHER_FUSIONCATCHER.out.versions)
        }

    emit:
        fusions  = ch_fusioncatcher_fusions     // channel [ meta, fusions ]
        versions = ch_versions                  // channel [ versions      ]
    }
