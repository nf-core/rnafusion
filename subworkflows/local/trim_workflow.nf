include { TRIMGALORE }                                  from '../../modules/nf-core/modules/trimgalore/main'

workflow TRIM_WORKFLOW {
    take:
        reads

    main:
        ch_versions = Channel.empty()
        // ch_dummy_file = file("$baseDir/assets/dummy_file_arriba.txt", checkIfExists: true)

        if (params.trim) {

            TRIMGALORE( reads )
            ch_versions = ch_versions.mix(TRIMGALORE.out.versions)
            ch_reads_out = TRIMGALORE.out.reads
        }
        else {
            ch_reads_out = reads
        }

    emit:
        reads         = ch_reads_out
        versions      = ch_versions.ifEmpty(null)
    }

