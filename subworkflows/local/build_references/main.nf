/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { GENCODE_DOWNLOAD }                from '../../../modules/local/gencode_download/main'
include { HGNC_DOWNLOAD }                   from '../../../modules/local/hgnc/main'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/
include { UCSC_GTFTOGENEPRED              } from '../../../modules/nf-core/ucsc/gtftogenepred/main'
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

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow BUILD_REFERENCES {

    take:
    // lists
    tools
    // strings
    fasta
    fai
    gtf
    genome_gencode_version
    genome
    hgnc_ref
    hgnc_date
    rrna_intervals
    refflat
    salmon_index
    starindex_ref
    arriba_ref_blacklist
    arriba_ref_cytobands
    arriba_ref_known_fusions
    arriba_ref_protein_domains
    fusioncatcher_ref
    starfusion_ref
    fusionreport_ref
    fusion_annot_lib
    species
    cosmic_username
    cosmic_passwd
    // paths
    pfam_file
    dfam_hmm
    dfam_h3f
    dfam_h3i
    dfam_h3m
    dfam_h3p
    annot_filter_url
    ctatsplicing_cancer_introns
    // boolean
    skip_vcf
    skip_qc
    no_cosmic

    main:
    def ch_versions = channel.empty()

    def ch_fasta = channel.empty()
    def ch_gtf   = channel.empty()
    if (!exists_not_empty(fasta) || !exists_not_empty(gtf)){
        GENCODE_DOWNLOAD(genome_gencode_version, genome)
        ch_fasta = GENCODE_DOWNLOAD.out.fasta.map { that -> [[id:that.Name], that] }
        ch_gtf = GENCODE_DOWNLOAD.out.gtf.map { that -> [[id:that.Name], that] }
    } else {
        ch_fasta = channel.fromPath(fasta).map { that -> [[id:that.Name], that] }
        ch_gtf = channel.fromPath(gtf).map { that -> [[id:that.Name], that] }
    }

    def ch_fai = channel.empty()
    if (!exists_not_empty(fai)){
        SAMTOOLS_FAIDX(ch_fasta.map { meta, fasta_ -> tuple(meta, fasta_, [])}, false)
        ch_fai = SAMTOOLS_FAIDX.out.fai
    } else {
        ch_fai = channel.fromPath(fai).map { that -> [[id:that.name.replaceFirst(/\.fai$/, '')], that] }
    }

    def ch_hgnc_date = channel.empty()
    def ch_hgnc_ref  = channel.empty()
    //TODO: unify as if(tools.contains("fusioninspector")) once nextflow bug fixed
    def run_fusioninspector = tools.contains("fusioninspector")
    if(run_fusioninspector && !skip_vcf) {
        if ((!exists_not_empty(hgnc_ref) || !exists_not_empty(hgnc_date)) && !skip_vcf){
            HGNC_DOWNLOAD( )
            ch_hgnc_ref = HGNC_DOWNLOAD.out.hgnc_ref.map { that -> [[id:that.name], that] }
            ch_hgnc_date = HGNC_DOWNLOAD.out.hgnc_date.map { that -> [[id:that.name], that] }
        } else {
            ch_hgnc_ref = channel.fromPath(hgnc_ref).map { that -> [[id:that.name], that] }
            ch_hgnc_date = channel.fromPath(hgnc_date).map { that -> [[id:that.name], that] }
        }
    }

    def ch_rrna_interval = channel.empty()
    if (!skip_qc) {
        if (!exists_not_empty(rrna_intervals)){
            GATK4_CREATESEQUENCEDICTIONARY(ch_fasta)

            BIOAWK(
                ch_gtf,
                [],
                false,
                "gff"
            )

            AGAT_CONVERTGFF2BED(BIOAWK.out.output)

            GATK4_BEDTOINTERVALLIST(AGAT_CONVERTGFF2BED.out.bed, GATK4_CREATESEQUENCEDICTIONARY.out.dict )

            ch_rrna_interval = GATK4_BEDTOINTERVALLIST.out.interval_list
        } else {
            ch_rrna_interval = channel.fromPath(rrna_intervals).map { that -> [[id:that.name], that] }
        }
    }

    def ch_refflat = channel.empty()
    if (!skip_qc) {
        if (!exists_not_empty(refflat)){
            UCSC_GTFTOGENEPRED(ch_gtf)
            ch_refflat = UCSC_GTFTOGENEPRED.out.refflat.map { meta, rf -> [[id: meta.id], rf] }
        } else {
            ch_refflat = channel.fromPath(refflat).map { that -> [[id:that.name], that] }
        }
    }

    def ch_salmon_index = channel.empty()
    if (tools.contains("salmon")) {
        if (!skip_qc) {
            if (!exists_not_empty(salmon_index)){
                GFFREAD(ch_gtf, ch_fasta.map{ _meta, fasta_ -> fasta_ })

                SALMON_INDEX(ch_fasta.map{ _meta, fasta_ -> fasta_ }, GFFREAD.out.gffread_fasta.map{ gffread_fasta -> gffread_fasta[1] })
                ch_salmon_index = SALMON_INDEX.out.index
            } else {
                ch_salmon_index = channel.fromPath(salmon_index)
            }
        }
    }

    def ch_starindex_ref = channel.empty()
    def star_index_tools = tools.intersect(["starfusion", "arriba", "ctatsplicing", "stringtie"])
    if (star_index_tools) {
        if (!exists_not_empty(starindex_ref)) {
            STAR_GENOMEGENERATE(ch_fasta, ch_gtf)
            ch_starindex_ref = STAR_GENOMEGENERATE.out.index
        } else {
            ch_starindex_ref = channel.fromPath(starindex_ref).map { that -> [[id:that.name], that] }
        }
    }

    def ch_arriba_ref_blacklist       = arriba_ref_blacklist ? channel.fromPath(arriba_ref_blacklist) : channel.empty()
    def ch_arriba_ref_cytobands       = arriba_ref_cytobands ? channel.fromPath(arriba_ref_cytobands) : channel.empty()
    def ch_arriba_ref_known_fusions   = arriba_ref_known_fusions ? channel.fromPath(arriba_ref_known_fusions) : channel.empty()
    def ch_arriba_ref_protein_domains = arriba_ref_protein_domains ? channel.fromPath(arriba_ref_protein_domains) : channel.empty()

    def ch_fusioncatcher_ref = fusioncatcher_ref ? channel.fromPath(fusioncatcher_ref).map { fusioncatcher_ref_ -> [[id:fusioncatcher_ref_.name], fusioncatcher_ref_] } : channel.empty()

    def ch_starfusion_ref = channel.empty()
    if (tools.intersect(["starfusion", "ctatsplicing", "fusioninspector"])) {
        if (!exists_not_empty(starfusion_ref)) {
            if(!fusion_annot_lib) {
                error("Expected --fusion_annot_lib to be specified when using StarFusion or any tools that depend on it")
            }

            if(pfam_file) {
                pfam_file = channel.value(pfam_file)
            } else {
                error("Expected `--pfam_version` to be specified when using StarFusion to automatically fill in Pfam database or specify `--pfam_file` for custom input")
            }

            if(dfam_hmm && dfam_h3p && dfam_h3m && dfam_h3i && dfam_h3f) {
                dfam_hmm = channel.value(dfam_hmm)
                dfam_h3f = channel.value(dfam_h3f)
                dfam_h3i = channel.value(dfam_h3i)
                dfam_h3m = channel.value(dfam_h3m)
                dfam_h3p = channel.value(dfam_h3p)
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

            STARFUSION_BUILD(ch_fasta, ch_gtf, fusion_annot_lib, species, pfam_file, dfam_urls_ch, annot_filter_url)
            ch_versions = ch_versions.mix(STARFUSION_BUILD.out.versions)
            if (tools.contains("ctatsplicing")) {
                CTATSPLICING_PREPGENOMELIB(
                    STARFUSION_BUILD.out.reference,
                    ctatsplicing_cancer_introns
                )
                ch_starfusion_ref = CTATSPLICING_PREPGENOMELIB.out.reference
            } else {
                ch_starfusion_ref = STARFUSION_BUILD.out.reference
            }
        }
        else {
            ch_starfusion_ref = channel.fromPath(starfusion_ref).map { starfusion_ref_ -> [[id:starfusion_ref_.name], starfusion_ref_] }
        }
    }

    def ch_fusionreport_ref = channel.empty()
    if (tools.contains("fusionreport")) {
        if (!exists_not_empty(fusionreport_ref)) {
            if (!no_cosmic && (!cosmic_username || !cosmic_passwd)) {
                error('COSMIC username and/or password missing, this is needed to download the fusionreport reference')
            }
            FUSIONREPORT_DOWNLOAD(channel.value([id:'fusionreport']))
            ch_fusionreport_ref = FUSIONREPORT_DOWNLOAD.out.fusionreport_ref
        } else {
            ch_fusionreport_ref = channel.fromPath(fusionreport_ref).map { that -> [[id:that.name], that] }
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
