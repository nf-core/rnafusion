/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { GENCODE_DOWNLOAD }                from '../../../modules/local/gencode_download/main'
include { HGNC_DOWNLOAD }                   from '../../../modules/local/hgnc/main'
include { GTF_TO_REFFLAT }                  from '../../../modules/local/uscs/custom_gtftogenepred/main'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/
include { CTATSPLICING_PREPGENOMELIB }      from '../../../modules/nf-core/ctatsplicing/prepgenomelib/main.nf'
include { BIOAWK                          } from '../../../modules/nf-core/bioawk/main'
include { AGAT_CONVERTGFF2BED             } from '../../../modules/nf-core/agat/convertgff2bed/main'
include { SAMTOOLS_FAIDX }                  from '../../../modules/nf-core/samtools/faidx/main'
include { STAR_GENOMEGENERATE }             from '../../../modules/nf-core/star/genomegenerate/main'
include { GATK4_CREATESEQUENCEDICTIONARY }  from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_BEDTOINTERVALLIST }         from '../../../modules/nf-core/gatk4/bedtointervallist/main'
include { SALMON_INDEX }                    from '../../../modules/nf-core/salmon/index/main'
include { FUSIONREPORT_DOWNLOAD }           from '../../../modules/nf-core/fusionreport/download/main'
include { STARFUSION_BUILD }                from '../../../modules/nf-core/starfusion/build/main'
include { GFFREAD }                         from '../../../modules/nf-core/gffread/main'
include { getFileSuffix } from '../../../modules/nf-core/cat/cat/main.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow BUILD_REFERENCES {

    take:
    tools // list of all the tools to create references for

    main:
    def ch_versions = Channel.empty()

    def ch_fasta = Channel.empty()
    def ch_gtf   = Channel.empty()
    if (!exists_not_empty(params.fasta) || !exists_not_empty(params.gtf)){
        GENCODE_DOWNLOAD(params.genome_gencode_version, params.genome)
        ch_versions = ch_versions.mix(GENCODE_DOWNLOAD.out.versions)
        ch_fasta = GENCODE_DOWNLOAD.out.fasta.map { that -> [[id:that.Name], that] }
        ch_gtf = GENCODE_DOWNLOAD.out.gtf.map { that -> [[id:that.Name], that] }
    } else {
        ch_fasta = Channel.fromPath(params.fasta).map { that -> [[id:that.Name], that] }
        ch_gtf = Channel.fromPath(params.gtf).map { that -> [[id:that.Name], that] }
    }

    def ch_fai = Channel.empty()
    if (!exists_not_empty(params.fai)){
        SAMTOOLS_FAIDX(ch_fasta, [[],[]], false)
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        ch_fai = SAMTOOLS_FAIDX.out.fai
    } else {
        ch_fai = Channel.fromPath(params.fai).map { that -> [[id:that.Name], that] }
    }

    def ch_hgnc_date = Channel.empty()
    def ch_hgnc_ref  = Channel.empty()
    //TODO: unify as if(tools.contains("fusioninspector")) once nextflow bug fixed
    def run_fusioninspector = tools.contains("fusioninspector")
    if(run_fusioninspector && !params.skip_vcf) {
        if ((!exists_not_empty(params.hgnc_ref) || !exists_not_empty(params.hgnc_date)) && !params.skip_vcf){
            HGNC_DOWNLOAD( )
            ch_versions = ch_versions.mix(HGNC_DOWNLOAD.out.versions)
            ch_hgnc_ref = HGNC_DOWNLOAD.out.hgnc_ref.map { that -> [[id:that.Name], that] }
            ch_hgnc_date = HGNC_DOWNLOAD.out.hgnc_date.map { that -> [[id:that.Name], that] }
        } else {
            ch_hgnc_ref = Channel.fromPath(params.hgnc_ref).map { that -> [[id:that.Name], that] }
            ch_hgnc_date = Channel.fromPath(params.hgnc_date).map { that -> [[id:that.Name], that] }
        }
    }

    def ch_rrna_interval = Channel.empty()
    if (!params.skip_qc) {
        if (!exists_not_empty(params.rrna_intervals)){
            GATK4_CREATESEQUENCEDICTIONARY(ch_fasta)
            ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)

            BIOAWK(ch_gtf)
            ch_versions = ch_versions.mix(BIOAWK.out.versions)

            AGAT_CONVERTGFF2BED(BIOAWK.out.output)
            ch_versions = ch_versions.mix(AGAT_CONVERTGFF2BED.out.versions)

            GATK4_BEDTOINTERVALLIST(AGAT_CONVERTGFF2BED.out.bed, GATK4_CREATESEQUENCEDICTIONARY.out.dict )
            ch_versions = ch_versions.mix(GATK4_BEDTOINTERVALLIST.out.versions)

            ch_rrna_interval = GATK4_BEDTOINTERVALLIST.out.interval_list
        } else {
            ch_rrna_interval = Channel.fromPath(params.rrna_intervals).map { that -> [[id:that.Name], that] }
        }
    }

    def ch_refflat = Channel.empty()
    if (!params.skip_qc) {
        if (!exists_not_empty(params.refflat)){
            GTF_TO_REFFLAT(ch_gtf)
            ch_versions = ch_versions.mix(GTF_TO_REFFLAT.out.versions)
            ch_refflat = GTF_TO_REFFLAT.out.refflat.map { that -> [[id:that.Name], that] }
        } else {
            ch_refflat = Channel.fromPath(params.refflat).map { that -> [[id:that.Name], that] }
        }
    }

    def ch_salmon_index = Channel.empty()
    if (tools.contains("salmon")) {
        if (!params.skip_qc) {
            if (!exists_not_empty(params.salmon_index)){
                GFFREAD(ch_gtf, ch_fasta.map{ it -> it[1] })
                ch_versions = ch_versions.mix(GFFREAD.out.versions)

                SALMON_INDEX(ch_fasta.map{ it -> it[1] }, GFFREAD.out.gffread_fasta.map{ it -> it[1] })
                ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)
                ch_salmon_index = SALMON_INDEX.out.index
            } else {
                ch_salmon_index = Channel.fromPath(params.salmon_index)
            }
        }
    }

    def ch_starindex_ref = Channel.empty()
    def star_index_tools = tools.intersect(["starfusion", "arriba", "ctatsplicing", "stringtie"])
    if (star_index_tools) {
        if (!exists_not_empty(params.starindex_ref)) {
            STAR_GENOMEGENERATE(ch_fasta, ch_gtf)
            ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
            ch_starindex_ref = STAR_GENOMEGENERATE.out.index
        } else {
            ch_starindex_ref = Channel.fromPath(params.starindex_ref).map { that -> [[id:that.Name], that] }
        }
    }

    def ch_arriba_ref_blacklist       = params.arriba_ref_blacklist ? Channel.fromPath(params.arriba_ref_blacklist) : Channel.empty()
    def ch_arriba_ref_cytobands       = params.arriba_ref_cytobands ? Channel.fromPath(params.arriba_ref_cytobands) : Channel.empty()
    def ch_arriba_ref_known_fusions   = params.arriba_ref_known_fusions ? Channel.fromPath(params.arriba_ref_known_fusions) : Channel.empty()
    def ch_arriba_ref_protein_domains = params.arriba_ref_protein_domains ? Channel.fromPath(params.arriba_ref_protein_domains) : Channel.empty()

    def ch_fusioncatcher_ref = params.fusioncatcher_ref ? Channel.fromPath(params.fusioncatcher_ref).map { it -> [[id:it.name], it] } : Channel.empty()

    def ch_starfusion_ref = Channel.empty()
    if (tools.intersect(["starfusion", "ctatsplicing", "fusioninspector"])) {
        if (exists_not_empty(params.starfusion_ref)) {
            ch_starfusion_ref = Channel.fromPath(params.starfusion_ref).map { it -> [[id:it.name], it] }
        }
        else {
            if(!params.fusion_annot_lib) {
                error("Expected --fusion_annot_lib to be specified when using StarFusion or any tools that depend on it")
            }

            if(params.pfam_file) {
                pfam_file = Channel.fromPath(params.pfam_file, checkIfExists: true)
            } else {
                error("Expected `--pfam_version` to be specified when using StarFusion to automatically fill in Pfam database or specify `--pfam_file` for custom input")
            }

            if(params.dfam_hmm && params.dfam_h3p && params.dfam_h3m && params.dfam_h3i && params.dfam_h3f) {
                dfam_hmm = Channel.fromPath(params.dfam_hmm, checkIfExists: true)
                dfam_h3f = Channel.fromPath(params.dfam_h3f, checkIfExists: true)
                dfam_h3i = Channel.fromPath(params.dfam_h3i, checkIfExists: true)
                dfam_h3m = Channel.fromPath(params.dfam_h3m, checkIfExists: true)
                dfam_h3p = Channel.fromPath(params.dfam_h3p, checkIfExists: true)
            } else {
                error("Expected `--dfam_version` and `--species` to be specified when using StarFusion to automatically fill in Dfam database or specify `--dfam_{hmm,h3f,h3i,h3m,h3p}` for custom input")
            }

            dfam_urls_ch = dfam_hmm
                .concat(
                    dfam_h3f,
                    dfam_h3i,
                    dfam_h3m,
                    dfam_h3p
                )
                .collect()

            STARFUSION_BUILD(ch_fasta, ch_gtf, params.fusion_annot_lib, params.species, pfam_file, dfam_urls_ch, params.annot_filter_url)
            ch_versions = ch_versions.mix(STARFUSION_BUILD.out.versions)
        }
        if (tools.contains("ctatsplicing") && !exists_not_empty(params.ctatsplicing_ref)) {
            CTATSPLICING_PREPGENOMELIB(
                ch_starfusion_ref,
                params.ctatsplicing_cancer_introns
            )
            ch_versions = ch_versions.mix(CTATSPLICING_PREPGENOMELIB.out.versions)
            ch_starfusion_ref = CTATSPLICING_PREPGENOMELIB.out.reference
        }
    }

    def ch_fusionreport_ref = Channel.empty()
    if (tools.contains("fusionreport")) {
        if (!exists_not_empty(params.fusionreport_ref)) {
            if (!params.no_cosmic && (!params.cosmic_username || !params.cosmic_passwd)) {
                error('COSMIC username and/or password missing, this is needed to download the fusionreport reference')
            }
            FUSIONREPORT_DOWNLOAD()
            ch_versions = ch_versions.mix(FUSIONREPORT_DOWNLOAD.out.versions)
            ch_fusionreport_ref = FUSIONREPORT_DOWNLOAD.out.fusionreport_ref
        } else {
            ch_fusionreport_ref = Channel.fromPath(params.fusionreport_ref).map { that -> [[id:that.Name], that] }
        }
    }

    emit:
    fasta                       = ch_fasta.collect()
    gtf                         = ch_gtf.collect()
    fai                         = ch_fai.collect()
    hgnc_ref                    = ch_hgnc_ref.collect()
    hgnc_date                   = ch_hgnc_date.collect()
    rrna_interval               = ch_rrna_interval.collect()
    refflat                     = ch_refflat.collect()
    salmon_index                = ch_salmon_index.collect()
    starindex_ref               = ch_starindex_ref.collect()
    arriba_ref_blacklist        = ch_arriba_ref_blacklist.collect()
    arriba_ref_cytobands        = ch_arriba_ref_cytobands.collect()
    arriba_ref_known_fusions    = ch_arriba_ref_known_fusions.collect()
    arriba_ref_protein_domains  = ch_arriba_ref_protein_domains.collect()
    fusioncatcher_ref           = ch_fusioncatcher_ref.collect()
    starfusion_ref              = ch_starfusion_ref.collect()
    fusionreport_ref            = ch_fusionreport_ref.collect()
    versions                    = ch_versions
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

//
// A function to test if a file exists and is not empty.
//   Input: A string that represents a file path
//   Output: A boolean
//
def exists_not_empty(path) {
    // Return false for invalid values
    if(!path) {
        return false
    }

    def path_to_check = file(path as String)
    // Return false if the path does not exist
    if(!path_to_check.exists()) {
        return false
    }

    // Don't check directories if the path is not local
    def is_local = path_to_check.getScheme() == "file"
    if(!is_local || !path_to_check.toFile().isDirectory()) {
        return !path_to_check.isEmpty()
    }

    // Get the first file in a directory and return whether it is empty or not
    def first_file = null
    path_to_check.toFile().eachFileRecurse(groovy.io.FileType.FILES) { file ->
        first_file = file
        return
    }
    return !first_file.toPath().isEmpty()
}
