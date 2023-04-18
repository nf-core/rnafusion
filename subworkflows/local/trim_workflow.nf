include { REFORMAT }                    from '../../modules/local/reformat/main'
include { FASTQC as FASTQC_FOR_TRIM }   from '../../modules/nf-core/fastqc/main'
include { FASTP }                       from '../../modules/nf-core/fastp/main'
include { FASTQC as FASTQC_FOR_FASTP }  from '../../modules/nf-core/fastqc/main'


workflow TRIM_WORKFLOW {
    take:
        reads

    main:
        ch_versions = Channel.empty()
        ch_reports = Channel.empty()

        if (params.trim) {

            REFORMAT( reads )
            ch_versions = ch_versions.mix(REFORMAT.out.versions)
            FASTQC_FOR_TRIM (REFORMAT.out.reads_out)
            ch_versions = ch_versions.mix(FASTQC_FOR_TRIM.out.versions)

            ch_reads_all = reads
            ch_reads_fusioncatcher = REFORMAT.out.reads_out
            ch_reports = FASTQC_FOR_TRIM.out.zip.collect{it[1]}.ifEmpty([])
        }
        else if (params.fastp_trim) {
            FASTP(reads, params.adapter_fasta, false, false)
            ch_versions = ch_versions.mix(FASTP.out.versions)

            FASTQC_FOR_FASTP(FASTP.out.reads)
            ch_versions = ch_versions.mix(FASTQC_FOR_FASTP.out.versions)

            ch_reads_all = FASTP.out.reads
            ch_reads_fusioncatcher = ch_reads_all
            ch_reports = ch_reports.mix(
                FASTQC_FOR_FASTP.out.zip.collect{it[1]}.ifEmpty([]),
                FASTP.out.json.collect{meta, json -> json},
                FASTP.out.html.collect{meta, html -> html}
            )
        }
        else {
            ch_reads_all = reads
            ch_reads_fusioncatcher = reads
        }

    emit:
        ch_reads_all
        ch_reads_fusioncatcher
        ch_reports
        versions      = ch_versions.ifEmpty(null)
    }

