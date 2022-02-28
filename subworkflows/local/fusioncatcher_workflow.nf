//
// Check input samplesheet and get read channels
//

include { FUSIONCATCHER }          from '../../modules/local/fusioncatcher/detect/main'
include { GET_PATH }               from '../../modules/local/getpath/main'


workflow FUSIONCATCHER_WORKFLOW {
    take:
        reads
        fast

    main:
        ch_versions = Channel.empty()
        ch_dummy_file = file("$baseDir/assets/dummy_file_fusioncatcher.txt", checkIfExists: true)

        if (params.fusioncatcher) {
            if (params.fusioncatcher_fusions){
                ch_fusioncatcher_fusions = params.fusioncatcher_fusions
            } else {
                FUSIONCATCHER (
                    INPUT_CHECK.out.reads,
                    params.fusioncatcher_ref
                )
                GET_PATH(FUSIONCATCHER.out.fusions)
                ch_fusioncatcher_fusions = GET_PATH.out.file
            }
        }
        else {
            ch_fusioncatcher_fusions = ch_dummy_file
        }

    emit:
        fusions  = ch_squid_fusions
        versions = ch_versions.ifEmpty(null)
    }

