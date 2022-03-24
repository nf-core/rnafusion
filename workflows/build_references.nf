/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { ENSEMBL_DOWNLOAD }                from '../modules/local/ensembl/main'
include { ARRIBA_DOWNLOAD }                 from '../modules/local/arriba/download/main'
include { ERICSCRIPT_DOWNLOAD }             from '../modules/local/ericscript/download/main'
include { FUSIONCATCHER_DOWNLOAD }          from '../modules/local/fusioncatcher/download/main'
include { FUSIONREPORT_DOWNLOAD }           from '../modules/local/fusionreport/download/main'
include { STARFUSION_DOWNLOAD }             from '../modules/local/starfusion/download/main'
include { GTF_TO_REFFLAT }                  from '../modules/local/uscs/custom_gtftogenepred/main'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

include { STAR_GENOMEGENERATE }             from '../modules/nf-core/modules/star/genomegenerate/main'
include { KALLISTO_INDEX as PIZZLY_INDEX }  from '../modules/nf-core/modules/kallisto/index/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow BUILD_REFERENCES {

    ENSEMBL_DOWNLOAD( params.ensembl_version )

    GTF_TO_REFFLAT(ENSEMBL_DOWNLOAD.out.chrgtf)

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

    if (params.fusionreport || params.all) {
        FUSIONREPORT_DOWNLOAD( params.cosmic_usr, params.cosmic_passwd )
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
