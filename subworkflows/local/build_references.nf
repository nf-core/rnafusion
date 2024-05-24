/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { ARRIBA_DOWNLOAD }                 from '../../modules/local/arriba/download/main'
include { ENSEMBL_DOWNLOAD }                from '../../modules/local/ensembl/main'
include { FUSIONCATCHER_DOWNLOAD }          from '../../modules/local/fusioncatcher/download/main'
include { FUSIONREPORT_DOWNLOAD }           from '../../modules/local/fusionreport/download/main'
include { HGNC_DOWNLOAD }                   from '../../modules/local/hgnc/main'
include { STARFUSION_BUILD }                from '../../modules/local/starfusion/build/main'
include { STARFUSION_DOWNLOAD }             from '../../modules/local/starfusion/download/main'
include { GTF_TO_REFFLAT }                  from '../../modules/local/uscs/custom_gtftogenepred/main'
include { RRNA_TRANSCRIPTS }                from '../../modules/local/rrnatranscripts/main'
include { CONVERT2BED }                     from '../../modules/local/convert2bed/main'
/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

include { SAMTOOLS_FAIDX }                  from '../../modules/nf-core/samtools/faidx/main'
include { STAR_GENOMEGENERATE }             from '../../modules/nf-core/star/genomegenerate/main'
include { GATK4_CREATESEQUENCEDICTIONARY }  from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_BEDTOINTERVALLIST }         from '../../modules/nf-core/gatk4/bedtointervallist/main'
include { SALMON_INDEX }                    from '../../modules/nf-core/salmon/index/main'
include { GFFREAD }                         from '../../modules/nf-core/gffread/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow BUILD_REFERENCES {


    main:

    ch_versions = Channel.empty()

    def fake_meta = [:]
    fake_meta.id = "Homo_sapiens.${params.genome}.${params.ensembl_version}"

    if (!file(params.fasta).exists() || file(params.fasta).isEmpty() ||
                    !file(params.chrgtf).exists() || file(params.chrgtf).isEmpty() ||
                    !file(params.gtf).exists() || file(params.gtf).isEmpty()){
        ENSEMBL_DOWNLOAD(params.ensembl_version, params.genome, fake_meta)}
    ch_fasta = {(!file(params.fasta).exists() || file(params.fasta).isEmpty()) ? ENSEMBL_DOWNLOAD.out.primary_assembly : Channel.fromPath(params.fasta).map { that -> [[id:that.Name], that] }.collect()}

    if (!file(params.hgnc_ref).exists() || file(params.hgnc_ref).isEmpty() ||
                    !file(params.hgnc_date).exists() || file(params.hgnc_date).isEmpty()){
        HGNC_DOWNLOAD( )}

    if (!file(params.fai).exists() || file(params.fai).isEmpty(){
        SAMTOOLS_FAIDX(ENSEMBL_DOWNLOAD.out.primary_assembly, [[],[]])}

    if (!file(params.rrna_intervals).exists() || file(params.rrna_intervals).isEmpty(){
        GATK4_CREATESEQUENCEDICTIONARY(ENSEMBL_DOWNLOAD.out.primary_assembly)
        RRNA_TRANSCRIPTS(ENSEMBL_DOWNLOAD.out.gtf)
        CONVERT2BED(RRNA_TRANSCRIPTS.out.rrna_gtf)
        GATK4_BEDTOINTERVALLIST(CONVERT2BED.out.bed, GATK4_CREATESEQUENCEDICTIONARY.out.dict)
        }

    if (!file(params.salmon_index).exists() || file(params.salmon_index).isEmpty(){ // add condition for qc, check that dirs can also be checked with isEmpty()
        GFFREAD(ENSEMBL_DOWNLOAD.out.gtf, ENSEMBL_DOWNLOAD.out.primary_assembly.map { meta, fasta -> [ fasta ] })
        SALMON_INDEX(ENSEMBL_DOWNLOAD.out.primary_assembly.map{ meta, fasta -> [ fasta ] }, GFFREAD.out.gffread_fasta.map{ meta, gffread_fasta -> [ gffread_fasta ] })
        }


    if ((params.starindex || params.all || params.starfusion || params.arriba) &&
            (!params.starindex_ref.exits() || params.starindex_ref.isEmpty())
            ) {
        STAR_GENOMEGENERATE( ENSEMBL_DOWNLOAD.out.primary_assembly, ENSEMBL_DOWNLOAD.out.gtf )
    }
    ch_starindex_ref = ...

    // if (params.arriba || params.all) {
    //     ARRIBA_DOWNLOAD()
    // }

    // if (params.fusioncatcher || params.all) {
    //     FUSIONCATCHER_DOWNLOAD()
    // }

    // if (params.starfusion || params.all) {
    //     if (params.starfusion_build){
    //         STARFUSION_BUILD( ENSEMBL_DOWNLOAD.out.primary_assembly, ENSEMBL_DOWNLOAD.out.gtf )
    //     } else {
    //         STARFUSION_DOWNLOAD()
    //     }
    // }

    // if (params.starfusion_build){
    //     GTF_TO_REFFLAT(ENSEMBL_DOWNLOAD.out.gtf)
    // } else {
    //     GTF_TO_REFFLAT(STARFUSION_DOWNLOAD.out.gtf)
    // }

    // if (params.fusionreport || params.all) {
    //     FUSIONREPORT_DOWNLOAD( params.cosmic_username, params.cosmic_passwd )
    // }

    emit:
    ch_fasta
    ch_chrgtf = {(!file(params.chrgtf).exists() || file(params.chrgtf).isEmpty()) ? ENSEMBL_DOWNLOAD.out.chrgtf : Channel.fromPath(params.chrgtf).map { that -> [[id:that.Name], that] }.collect()}
    ch_gtf = {(!file(params.gtf).exists() || file(params.gtf).isEmpty()) ? ENSEMBL_DOWNLOAD.out.gtf : Channel.fromPath(params.gtf).map { that -> [[id:that.Name], that] }.collect()}
    ch_hgnc_ref = Channel.fromPath(params.hgnc_ref).map { it -> [[id:it.Name], it] }.collect()
    ch_hgnc_date = Channel.fromPath(params.hgnc_date).map { it -> [[id:it.Name], it] }.collect()
    ch_fai = Channel.fromPath(params.fai).map { it -> [[id:it.Name], it] }.collect()
    ch_rrna_interval = params.starfusion_build ?  Channel.fromPath(params.rrna_intervals).map { it -> [[id:it.Name], it] }.collect() : Channel.fromPath("${params.ensembl_ref}/ref_annot.interval_list").map { it -> [[id:it.Name], it] }.collect()
    ch_salmon_index = Channel.fromPath(params.salmon_index).map { it -> [[id:it.Name], it] }.collect()


}

/*
========================================================================================
    THE END
========================================================================================
*/
