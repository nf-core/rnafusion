include { FUSIONCATCHER }                       from '../../modules/local/fusioncatcher/detect/main'


workflow FUSIONCATCHER_WORKFLOW {
    take:
        reads

    main:
        ch_versions = Channel.empty()
        ch_dummy_file = file("$baseDir/assets/dummy_file_fusioncatcher.txt", checkIfExists: true)

        if ((params.fusioncatcher || params.all) && !params.fusioninspector_only) {
            if (params.fusioncatcher_fusions){
                ch_fusioncatcher_fusions = reads.combine(Channel.value(file(params.fusioncatcher_fusions, checkIfExists:true)))
                                            .map { meta, reads, fusions -> [ meta, fusions ] }
            } else {
                FUSIONCATCHER (
                    reads,
                    params.fusioncatcher_ref
                )
                ch_fusioncatcher_fusions = FUSIONCATCHER.out.fusions
                ch_versions = ch_versions.mix(FUSIONCATCHER.out.versions)
            }
        }
        else {
            ch_fusioncatcher_fusions = reads.combine(Channel.value(file(ch_dummy_file, checkIfExists:true)))
                                        .map { meta, reads, fusions -> [ meta, fusions ] }
        }

    emit:
        fusions  = ch_fusioncatcher_fusions
        versions = ch_versions
    }

