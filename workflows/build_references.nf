/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { ENSEMBL_DOWNLOAD }                from '../modules/local/ensembl/main'
include { FUSIONCATCHER_DOWNLOAD }          from '../modules/local/fusioncatcher/download/main'
include { FUSIONREPORT_DOWNLOAD }           from '../modules/local/fusionreport/download/main'
include { HGNC_DOWNLOAD }                   from '../modules/local/hgnc/main'
include { STARFUSION_BUILD }                from '../modules/local/starfusion/build/main'
include { STARFUSION_DOWNLOAD }             from '../modules/local/starfusion/download/main'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

include { ARRIBA_DOWNLOAD }                 from '../modules/nf-core/arriba/download/main'
include { BEDOPS_CONVERT2BED }              from '../modules/nf-core/bedops/convert2bed/main'
include { GATK4_BEDTOINTERVALLIST }         from '../modules/nf-core/gatk4/bedtointervallist/main'
include { GATK4_CREATESEQUENCEDICTIONARY }  from '../modules/nf-core/gatk4/createsequencedictionary/main'
include { RRNATRANSCRIPTS }                 from '../modules/nf-core/rrnatranscripts/main'
include { SAMTOOLS_FAIDX }                  from '../modules/nf-core/samtools/faidx/main'
include { STAR_GENOMEGENERATE }             from '../modules/nf-core/star/genomegenerate/main'
include { UCSC_GTFTOGENEPRED }              from '../modules/nf-core/ucsc/gtftogenepred/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow BUILD_REFERENCES {

    def fake_meta = [:]
    fake_meta.id = "Homo_sapiens.${params.genome}.${params.ensembl_version}"
    ENSEMBL_DOWNLOAD( params.ensembl_version, params.genome, fake_meta )
    HGNC_DOWNLOAD()


    SAMTOOLS_FAIDX(ENSEMBL_DOWNLOAD.out.fasta, [[],[]])
    GATK4_CREATESEQUENCEDICTIONARY(ENSEMBL_DOWNLOAD.out.fasta)

    RRNATRANSCRIPTS(ENSEMBL_DOWNLOAD.out.gtf.map{ meta, gtf -> [ gtf ] })
    BEDOPS_CONVERT2BED(RRNATRANSCRIPTS.out.rrna_gtf.map{ it -> [[id:it.Name], it] })
    GATK4_BEDTOINTERVALLIST(BEDOPS_CONVERT2BED.out.bed, GATK4_CREATESEQUENCEDICTIONARY.out.dict)


    if (params.starindex || params.all || params.starfusion || params.arriba) {
        STAR_GENOMEGENERATE( ENSEMBL_DOWNLOAD.out.fasta, ENSEMBL_DOWNLOAD.out.gtf )
    }

    if (params.arriba || params.all) {
        ARRIBA_DOWNLOAD()
    }

    if (params.fusioncatcher || params.all) {
        FUSIONCATCHER_DOWNLOAD()
    }

    if (params.starfusion || params.all) {
        if (params.starfusion_build){
            STARFUSION_BUILD( ENSEMBL_DOWNLOAD.out.fasta, ENSEMBL_DOWNLOAD.out.chrgtf )
        } else {
            STARFUSION_DOWNLOAD()
        }
    }

    if (params.starfusion_build){
        UCSC_GTFTOGENEPRED(ENSEMBL_DOWNLOAD.out.chrgtf)
    } else {
        UCSC_GTFTOGENEPRED(STARFUSION_DOWNLOAD.out.chrgtf)
    }

    if (params.fusionreport || params.all) {
        FUSIONREPORT_DOWNLOAD( params.cosmic_username, params.cosmic_passwd )
    }

}

/*
========================================================================================
    THE END
========================================================================================
*/
