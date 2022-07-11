include { REFORMAT }                    from '../../modules/local/reformat/main'
include { FASTQC as FASTQC_FOR_TRIM }   from '../modules/nf-core/modules/fastqc/main'

workflow TRIM_WORKFLOW {
    take:
        reads

    main:
        ch_versions = Channel.empty()

        if (params.trim) {

            REFORMAT( reads )
            ch_versions = ch_versions.mix(REFORMAT.out.versions)
            ch_reads_out = REFORMAT.out.reads_out

            FASTQC_FOR_TRIM (ch_reads_out)
        }
        else {
            ch_reads_out = reads
        }

    emit:
        reads         = ch_reads_out
        versions      = ch_versions.ifEmpty(null)
    }

