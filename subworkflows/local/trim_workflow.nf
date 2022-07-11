include { REFORMAT }                 from '../../modules/local/reformat/main'

workflow TRIM_WORKFLOW {
    take:
        reads

    main:
        ch_versions = Channel.empty()
        // ch_dummy_file = file("$baseDir/assets/dummy_file_arriba.txt", checkIfExists: true)

        if (params.trim) {

            REFORMAT( reads )
            ch_versions = ch_versions.mix(REFORMAT.out.versions)
            ch_reads_out = REFORMAT.out.reads_out
        }
        else {
            ch_reads_out = reads
        }

    emit:
        reads         = ch_reads_out
        versions      = ch_versions.ifEmpty(null)
    }

