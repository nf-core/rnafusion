include { REFORMAT }                    from '../../modules/local/reformat/main'
include { FASTQC as FASTQC_FOR_TRIM }   from '../../modules/nf-core/fastqc/main'
include { FASTP }                       from '../../modules/nf-core/fastp/main'
include { FASTQC as FASTQC_FOR_FASTP }  from '../../modules/nf-core/fastqc/main'


workflow TRIM_WORKFLOW {
    take:
        reads

    main:
        ch_versions = Channel.empty()

        if (params.trim) {

            REFORMAT( reads )
            ch_versions = ch_versions.mix(REFORMAT.out.versions)
            FASTQC_FOR_TRIM (REFORMAT.out.reads_out)
            ch_versions = ch_versions.mix(FASTQC_FOR_TRIM.out.versions)

            ch_reads_all = reads
            ch_reads_fusioncatcher = REFORMAT.out.reads_out
        }
        else if (params.fastp_trim) {
            FASTP(reads, params.adapter_fasta, false, false)
            ch_versions = ch_versions.mix(FASTP.out.versions)

            FASTQC_FOR_FASTP(FASTP.out.reads)
            ch_versions = ch_versions.mix(FASTQC_FOR_FASTP.out.versions)

            ch_reads_all = FASTP.out.reads
            ch_reads_fusioncatcher = ch_reads_all

        }
        else {
            ch_reads_all = reads
            ch_reads_fusioncatcher = reads
        }

    emit:
        ch_reads_all
        ch_reads_fusioncatcher
        versions      = ch_versions.ifEmpty(null)
    }

