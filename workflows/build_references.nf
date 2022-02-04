/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { ENSEMBL_DOWNLOAD }                from '../modules/local/ensembl/main'
include { ARRIBA_DOWNLOAD }                 from '../modules/local/arriba/download/main'
include { STAR_GENOMEGENERATE }             from '../modules/local/star/genomegenerate/main'
include { ERICSCRIPT_DOWNLOAD }             from '../modules/local/ericscript/download/main'
include { FUSIONCATCHER_DOWNLOAD }          from '../modules/local/fusioncatcher/download/main'
include { KALLISTO_INDEX as PIZZLY_INDEX }  from '../modules/nf-core/modules/kallisto/index/main'
include { STARFUSION_DOWNLOAD }             from '../modules/local/starfusion/download/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow BUILD_REFERENCES {

    ENSEMBL_DOWNLOAD( params.ensembl_version )

    if (params.starindex || params.all) {
        STAR_GENOMEGENERATE( ENSEMBL_DOWNLOAD.out.fasta, ENSEMBL_DOWNLOAD.out.gtf )
    }

    if (params.arriba || params.all) {
        ARRIBA_DOWNLOAD()
    }

    if (params.ericscript || params.all) {
        ERICSCRIPT_DOWNLOAD()
    }

    if (params.fusioncatcher || params.all) {
        FUSIONCATCHER_DOWNLOAD()
    }

    if (params.pizzly || params.all) {
        PIZZLY_INDEX( ENSEMBL_DOWNLOAD.out.transcript )
    }

    if (params.starfusion || params.all) {
        STARFUSION_DOWNLOAD( ENSEMBL_DOWNLOAD.out.fasta, ENSEMBL_DOWNLOAD.out.chrgtf )
    }

}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
    NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
