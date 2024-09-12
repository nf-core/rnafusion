include { FASTP }                       from '../../modules/nf-core/fastp/main'
include { FASTQC as FASTQC_FOR_FASTP }  from '../../modules/nf-core/fastqc/main'


workflow TRIM_WORKFLOW {
    take:
        reads

    main:
        ch_versions = Channel.empty()

        if (params.fastp_trim) {
            FASTP(reads, params.adapter_fasta, false, false)
            ch_versions = ch_versions.mix(FASTP.out.versions)

            FASTQC_FOR_FASTP(FASTP.out.reads)
            ch_versions = ch_versions.mix(FASTQC_FOR_FASTP.out.versions)

            ch_reads = FASTP.out.reads
            ch_fastp_html = FASTP.out.html
            ch_fastp_json = FASTP.out.json
            ch_fastqc_trimmed = FASTQC_FOR_FASTP.out.zip

        }
        else {
            ch_reads = reads
        }

    emit:
        trimmed_reads = ch_reads
        fastp_html = ch_fastp_html.ifEmpty([])
        fastp_json = ch_fastp_json.ifEmpty([])
        fastqc_trimmed = ch_fastqc_trimmed.ifEmpty([])
        versions = ch_versions
    }

