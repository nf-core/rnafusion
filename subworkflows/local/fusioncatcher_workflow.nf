include { FUSIONCATCHER }                       from '../../modules/local/fusioncatcher/detect/main'


workflow FUSIONCATCHER_WORKFLOW {
    take:
        reads

    main:
        ch_versions = Channel.empty()
        ch_dummy_file = file("$baseDir/assets/dummy_file_fusioncatcher.txt", checkIfExists: true)

        if ((params.fusioncatcher || params.all) && !params.fusioninspector_only) {
            if (params.fusioncatcher_fusions){
                ch_fusioncatcher_fusions = reads.merge(Channel.fromPath(params.fusioncatcher_fusions, checkIfExists:true))
                
            } else {
                FUSIONCATCHER (
                    reads,
                    params.fusioncatcher_ref
                )
                ch_fusioncatcher_fusions = FUSIONCATCHER.out.fusions
            }
        }
        else {
            ch_fusioncatcher_fusions = reads.merge(ch_dummy_file)
        }

    emit:
        fusions  = ch_fusioncatcher_fusions
        versions = ch_versions.ifEmpty(null)
    }

