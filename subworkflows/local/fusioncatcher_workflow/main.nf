include { FASTP as FASTP_FOR_FUSIONCATCHER } from '../../../modules/nf-core/fastp/main'
include { FUSIONCATCHER_FUSIONCATCHER }      from '../../../modules/nf-core/fusioncatcher/fusioncatcher/main'

// TODO: Remove fusioncatcher_fusions as parameter.
// TODO: remove dummy file. Work with channel.empty()
// TODO: if the files were already produced and the user want to skip the module because of this, they should be taken them from the sample sheet

workflow FUSIONCATCHER_WORKFLOW {
    take:
        reads                   // channel [ meta, [ fastqs ] ]
        fusioncatcher_trimming  // boolean
        adapter_fasta           // path
        fusioncatcher_ref       // channel [ meta, path       ]
        fusioncatcher_fusions   // path, string

    main:
        ch_versions   = channel.empty()

        if (fusioncatcher_fusions){

            ch_fusioncatcher_fusions = reads.combine(channel.value(file(fusioncatcher_fusions, checkIfExists:true)))
                                        .map { meta, _reads, fusions -> [ meta, fusions ] }
        } else {
            if (fusioncatcher_trimming) {
                reads_with_adapters = reads
                    .map { meta, reads_files ->
                        [ meta, reads_files, adapter_fasta ? file(adapter_fasta, checkIfExists: true) : [] ]
                    }

                FASTP_FOR_FUSIONCATCHER(
                    reads_with_adapters,
                    false, // discard_trimmed_pass
                    false, // save_trimmed_fail
                    false  // save_merged
                )

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
