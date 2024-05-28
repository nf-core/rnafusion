/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { ARRIBA_DOWNLOAD }                 from '../../modules/local/arriba/download/main'
include { ENSEMBL_DOWNLOAD }                from '../../modules/local/ensembl/main'
include { FUSIONCATCHER_DOWNLOAD }          from '../../modules/local/fusioncatcher/download/main'
include { FUSIONCATCHER_BUILD }             from '../../modules/local/fusioncatcher/build/main'
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


    if (!file(params.fasta).exists() || file(params.fasta).isEmpty() ||
            !file(params.chrgtf).exists() || file(params.chrgtf).isEmpty() ||
            !file(params.gtf).exists() || file(params.gtf).isEmpty()){
        fake_meta = [:]
        fake_meta.id = "Homo_sapiens.${params.genome}.${params.ensembl_version}"
        ENSEMBL_DOWNLOAD(params.ensembl_version, params.genome, fake_meta)}
    ch_fasta = Channel.fromPath(params.fasta).map { that -> [[id:that.Name], that] }
    ch_chrgtf = Channel.fromPath(params.chrgtf).map { that -> [[id:that.Name], that] }
    ch_gtf = Channel.fromPath(params.gtf).map { that -> [[id:that.Name], that] }

    if (!file(params.hgnc_ref).exists() || file(params.hgnc_ref).isEmpty() ||
                !file(params.hgnc_date).exists() || file(params.hgnc_date).isEmpty()){
        HGNC_DOWNLOAD( )}
    ch_hgnc_ref = Channel.fromPath(params.hgnc_ref).map { that -> [[id:that.Name], that] }
    ch_hgnc_date = Channel.fromPath(params.hgnc_date).map { that -> [[id:that.Name], that] }


    if (!file(params.fai).exists() || file(params.fai).isEmpty()){
        SAMTOOLS_FAIDX(ch_fasta, [[],[]])}
    ch_fai = Channel.fromPath(params.fai).map { that -> [[id:that.Name], that] }


    if (!file(params.rrna_intervals).exists() || file(params.rrna_intervals).isEmpty()){
        GATK4_CREATESEQUENCEDICTIONARY(ch_fasta)
        RRNA_TRANSCRIPTS(ch_gtf)
        CONVERT2BED(RRNA_TRANSCRIPTS.out.rrna_gtf)
        GATK4_BEDTOINTERVALLIST(CONVERT2BED.out.bed, GATK4_CREATESEQUENCEDICTIONARY.out.dict)
        }
    ch_rrna_interval = Channel.fromPath(params.rrna_intervals).map { that -> [[id:that.Name], that] }


    if (!file(params.refflat).exists() || file(params.refflat).isEmpty()){
        GTF_TO_REFFLAT(ch_gtf)}
    ch_refflat = Channel.fromPath(params.refflat)


    if (!file(params.salmon_index).exists() || file(params.salmon_index).isEmpty() ||
        !file(params.salmon_index_stub_check).exists() || file(params.salmon_index_stub_check).isEmpty()){ // add condition for qc
        GFFREAD(ch_gtf, ch_fasta.map{ meta, fasta -> [ fasta ] })
        SALMON_INDEX(ch_fasta.map{ meta, fasta -> [ fasta ] }, GFFREAD.out.gffread_fasta.map{ meta, gffread_fasta -> [ gffread_fasta ] })
        }
    ch_salmon_index = Channel.fromPath(params.salmon_index)


    if ((params.starindex || params.all || params.starfusion || params.arriba) &&
            (!file(params.starindex_ref).exists() || file(params.starindex_ref).isEmpty() ||
            !file(params.starindex_ref_stub_check).exists() || file(params.starindex_ref_stub_check).isEmpty() )) {
        STAR_GENOMEGENERATE(ch_fasta, ch_gtf)
    }
    ch_starindex_ref = Channel.fromPath(params.starindex_ref).map { that -> [[id:that.Name], that] }


    if ((params.arriba || params.all) &&
            (!file(params.arriba_ref_blacklist).exists() || file(params.arriba_ref_blacklist).isEmpty() ||
            !file(params.arriba_ref_known_fusions).exists() || file(params.arriba_ref_known_fusions).isEmpty() ||
            !file(params.arriba_ref_protein_domains).exists() || file(params.arriba_ref_protein_domains).isEmpty())) {
        ARRIBA_DOWNLOAD()}
    ch_arriba_ref_blacklist = Channel.fromPath(params.arriba_ref_blacklist).map { that -> [[id:that.Name], that] }
    ch_arriba_ref_known_fusions = Channel.fromPath(params.arriba_ref_known_fusions).map { that -> [[id:that.Name], that] }
    ch_arriba_ref_protein_domains = Channel.fromPath(params.arriba_ref_protein_domains).map { that -> [[id:that.Name], that] }


    if ((params.fusioncatcher || params.all) &&
            (!file(params.fusioncatcher_ref).exists() || file(params.fusioncatcher_ref).isEmpty() ||
            !file(params.fusioncatcher_ref_stub_check).exists() || file(params.fusioncatcher_ref_stub_check).isEmpty() )) {
        if (params.download_refs) {
            FUSIONCATCHER_DOWNLOAD(params.ensembl_version)}
        else {
            FUSIONCATCHER_BUILD(params.ensembl_version)}
    }
    ch_fusioncatcher_ref = Channel.fromPath(params.fusioncatcher_ref).map { that -> [[id:that.Name], that] }


    if ((params.starfusion || params.all) &&
            (!file(params.starfusion_ref).exists() || file(params.starfusion_ref).isEmpty() ||
            !file(params.starfusion_ref_stub_check).exists() || file(params.starfusion_ref_stub_check).isEmpty() )) {
        if (params.download_refs) {
            STARFUSION_DOWNLOAD( ch_fasta, ch_gtf )}
        else {
            STARFUSION_BUILD( ch_fasta, ch_gtf )}
    }
    ch_starfusion_ref = Channel.fromPath(params.starfusion_ref).map { that -> [[id:that.Name], that] }


    if ((params.fusionreport || params.all) &&
            (!file(params.fusionreport_ref).exists() || file(params.fusionreport_ref).isEmpty() ||
            !file(params.fusionreport_ref_stub_check).exists() || file(params.fusionreport_ref_stub_check).isEmpty())) {
        if (!params.cosmic_username || !params.cosmic_passwd) { exit 1, 'COSMIC username and/or password missing' }
        FUSIONREPORT_DOWNLOAD( params.cosmic_username, params.cosmic_passwd )}
    ch_fusionreport_ref = Channel.fromPath(params.fusionreport_ref).map { that -> [[id:that.Name], that] }

    emit:
    ch_fasta
    ch_chrgtf
    ch_gtf
    ch_hgnc_ref
    ch_hgnc_date
    ch_fai
    ch_rrna_interval
    ch_refflat
    ch_salmon_index
    ch_starindex_ref
    ch_arriba_ref_blacklist
    ch_arriba_ref_known_fusions
    ch_arriba_ref_protein_domains
    ch_fusioncatcher_ref
    ch_starfusion_ref
    ch_fusionreport_ref
}

/*
========================================================================================
    THE END
========================================================================================
*/
