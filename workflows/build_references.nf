/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { ARRIBA_DOWNLOAD }                 from '../modules/local/arriba/download/main'
include { ENSEMBL_DOWNLOAD }                from '../modules/local/ensembl/main'
include { FUSIONCATCHER_DOWNLOAD }          from '../modules/local/fusioncatcher/download/main'
include { FUSIONREPORT_DOWNLOAD }           from '../modules/local/fusionreport/download/main'
include { HGNC_DOWNLOAD }                   from '../modules/local/hgnc/main'
include { STARFUSION_BUILD }                from '../modules/local/starfusion/build/main'
include { STARFUSION_DOWNLOAD }             from '../modules/local/starfusion/download/main'
include { GTF_TO_REFFLAT }                  from '../modules/local/uscs/custom_gtftogenepred/main'
include { RRNA_TRANSCRIPTS }                from '../modules/local/rrnatranscripts/main'
include { CONVERT2BED }                     from '../modules/local/convert2bed/main'
/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

include { SAMTOOLS_FAIDX }                  from '../modules/nf-core/samtools/faidx/main'
include { STAR_GENOMEGENERATE }             from '../modules/nf-core/star/genomegenerate/main'
include { GATK4_CREATESEQUENCEDICTIONARY }  from '../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_BEDTOINTERVALLIST }         from '../modules/nf-core/gatk4/bedtointervallist/main'
include { SALMON_INDEX }                    from '../modules/nf-core/salmon/index/main'
include { GFFREAD                                      } from '../modules/local/gffread'
include { GUNZIP } from '../modules/nf-core/gunzip/main'
/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow BUILD_REFERENCES {

    def fake_meta = [:]
    fake_meta.id = "Homo_sapiens.${params.genome}.${params.ensembl_version}"
    ENSEMBL_DOWNLOAD( params.ensembl_version, params.genome, fake_meta )
    // HGNC_DOWNLOAD( )
    GUNZIP(ENSEMBL_DOWNLOAD.out.primary_assembly)


    SAMTOOLS_FAIDX(GUNZIP.out.gunzip, [[],[]])
    // GATK4_CREATESEQUENCEDICTIONARY(ENSEMBL_DOWNLOAD.out.fasta)

    // RRNA_TRANSCRIPTS(ENSEMBL_DOWNLOAD.out.gtf)
    // CONVERT2BED(RRNA_TRANSCRIPTS.out.rrna_gtf)

    // GATK4_BEDTOINTERVALLIST(CONVERT2BED.out.bed, GATK4_CREATESEQUENCEDICTIONARY.out.dict)

    ch_fasta_fai =(GUNZIP.out.gunzip).join(SAMTOOLS_FAIDX.out.fai)
    GFFREAD(ENSEMBL_DOWNLOAD.out.gtf, ch_fasta_fai)
    SALMON_INDEX(GUNZIP.out.gunzip.map { meta, fasta -> [ fasta ] }, GFFREAD.out.tr_fasta)

    // RSEM_PREPAREREFERENCE(ENSEMBL_DOWNLOAD.out.fasta.map { meta, fasta -> [ fasta ] }, ENSEMBL_DOWNLOAD.out.gtf.map { meta, gtf -> [ gtf ] })


    // if (params.starindex || params.all || params.starfusion || params.arriba) {
    //     STAR_GENOMEGENERATE( ENSEMBL_DOWNLOAD.out.fasta, ENSEMBL_DOWNLOAD.out.gtf )
    // }

    // if (params.arriba || params.all) {
    //     ARRIBA_DOWNLOAD()
    // }

    // if (params.fusioncatcher || params.all) {
    //     FUSIONCATCHER_DOWNLOAD()
    // }

    // if (params.starfusion || params.all) {
    //     if (params.starfusion_build){
    //         STARFUSION_BUILD( ENSEMBL_DOWNLOAD.out.fasta, ENSEMBL_DOWNLOAD.out.chrgtf )
    //     } else {
    //         STARFUSION_DOWNLOAD()
    //     }
    // }

    // if (params.starfusion_build){
    //     GTF_TO_REFFLAT(ENSEMBL_DOWNLOAD.out.chrgtf)
    // } else {
    //     GTF_TO_REFFLAT(STARFUSION_DOWNLOAD.out.chrgtf)
    // }

    // if (params.fusionreport || params.all) {
    //     FUSIONREPORT_DOWNLOAD( params.cosmic_username, params.cosmic_passwd )
    // }

}

/*
========================================================================================
    THE END
========================================================================================
*/
