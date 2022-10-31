/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { ENSEMBL_DOWNLOAD }                from '../modules/local/ensembl/main'
include { FUSIONCATCHER_DOWNLOAD }          from '../modules/local/fusioncatcher/download/main'
include { FUSIONREPORT_DOWNLOAD }           from '../modules/local/fusionreport/download/main'
include { STARFUSION_BUILD }                from '../modules/local/starfusion/build/main'
include { STARFUSION_DOWNLOAD }             from '../modules/local/starfusion/download/main'
include { GTF_TO_REFFLAT }                  from '../modules/local/uscs/custom_gtftogenepred/main'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

include { SAMTOOLS_FAIDX }                  from '../modules/nf-core/samtools/faidx/main'
include { STAR_GENOMEGENERATE }             from '../modules/nf-core/star/genomegenerate/main'
include { KALLISTO_INDEX as PIZZLY_INDEX }  from '../modules/nf-core/kallisto/index/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow BUILD_REFERENCES {

    ENSEMBL_DOWNLOAD( params.ensembl_version )
    // ch_meta_empty = Channel.value([[:])

    // ch_fasta_path = Channel.from(ENSEMBL_DOWNLOAD.out.fasta)
    // ch_meta_empty = Channel.value([:]).view()
    ch_fasta_w_meta = ENSEMBL_DOWNLOAD.out.fasta.map{ it -> [[id:it[0].baseName], it] }
    SAMTOOLS_FAIDX(ch_fasta_w_meta)


    if (params.starindex || params.all || params.starfusion || params.arriba || params.squid ) {
        STAR_GENOMEGENERATE( ENSEMBL_DOWNLOAD.out.fasta, ENSEMBL_DOWNLOAD.out.gtf )
    }

    if (params.fusioncatcher || params.all) {
        FUSIONCATCHER_DOWNLOAD()
    }

    if (params.pizzly || params.all) {
        PIZZLY_INDEX( ENSEMBL_DOWNLOAD.out.transcript )
    }

    if (params.starfusion || params.all) {
        if (params.starfusion_build){
            STARFUSION_BUILD( ENSEMBL_DOWNLOAD.out.fasta, ENSEMBL_DOWNLOAD.out.chrgtf )
        } else {
            STARFUSION_DOWNLOAD()
        }
    }

    if (params.starfusion_build){
        GTF_TO_REFFLAT(ENSEMBL_DOWNLOAD.out.chrgtf)
    } else {
        GTF_TO_REFFLAT(STARFUSION_DOWNLOAD.out.chrgtf)
    }

    if (params.fusionreport || params.all) {
        FUSIONREPORT_DOWNLOAD( params.cosmic_username, params.cosmic_passwd )
    }

}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
