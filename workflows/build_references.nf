/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include { ENSEMBL_DOWNLOAD }                from '../modules/local/ensembl/main'                    addParams( options: modules['ensembl_download'] )
include { ARRIBA_DOWNLOAD }                 from '../modules/local/arriba/download/main'            addParams( options: modules['arriba_download'] )
include { STAR_GENOMEGENERATE }             from '../modules/local/star/genomegenerate/main'        addParams( options: modules['star_index'] )
include { ERICSCRIPT_DOWNLOAD }             from '../modules/local/ericscript/download/main'        addParams( options: modules['ericscript_download'] )
include { FUSIONCATCHER_DOWNLOAD }          from '../modules/local/fusioncatcher/download/main'     addParams( options: modules['fusioncatcher_download'] )
include { KALLISTO_INDEX as PIZZLY_INDEX }  from '../modules/nf-core/modules/kallisto/index/main'   addParams( options: modules['pizzly_download'] )
include { STARFUSION_DOWNLOAD }             from '../modules/local/starfusion/download/main'        addParams( options: modules['starfusion_download'] )

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
        STARFUSION_DOWNLOAD( ENSEMBL_DOWNLOAD.out.fasta, ENSEMBL_DOWNLOAD.out.gtf )
    }

}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
    // NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/