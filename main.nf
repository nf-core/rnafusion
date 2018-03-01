#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
               N G I - R N A S E Q    F U S I O N D E T E C T
========================================================================================
 New RNA-Seq Best Practice Analysis Pipeline. Started May 2017.
 #### Homepage / Documentation
 https://github.com/SciLifeLab/NGI-RNAfusion
 #### Authors
 Rickard Hammarén @Hammarn  <rickard.hammaren@scilifelab.se>
 Philip Ewels @ewels <phil.ewels@scilifelab.se>
*/

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = '0.1'

// Configurable variables - same as NGI-RNAseq for now
params.project = false
params.reads = "data/*{1,2}.fastq.gz"
params.email = false
params.star_fusion = true
params.inspector = false
params.fusioncatcher = true
params.sensitivity = 'sensitive'
params.clusterOptions = false
params.outdir = './results'
params.fc_extra_options = ''

params.singleEnd = false
Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 ) 
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}" }
    .into { read_files_star_fusion; fusion_inspector_reads; fusioncatcher_reads}


// Validate inputs
if( params.star_fusion_reference && params.star_fusion ){
    Channel.fromPath(params.star_fusion_reference)
        .ifEmpty { exit 1, "STAR-fusion reference not found: ${params.star_fusion_reference}" }
        .into { star_fusion_reference, star_fusion_reference_fusioninspector }
    }
    if( params.fusioncatcher_data_dir && params.fusioncatcher ){
        fusioncatcher_data_dir = Channel
        .fromPath(params.fusioncatcher_data_dir)
        .ifEmpty { exit 1, "FusionCatcher data directory not found: ${params.fusioncatcher_data_dir}" }
}

/*
 * STAR-Fusion
 */
process star_fusion{
    tag "$name"
    publishDir "${params.outdir}/Star_fusion",  mode: 'copy'
    
    input:
    set val (name), file(reads) from read_files_star_fusion
    file star_fusion_reference from star_fusion_reference.collect()
    output:
    file '*final.abridged*' into star_fusion_abridged
    file '*star_fusion.fusion_candidates.final.abridged.FFPM' into fusion_candidates,fusion_candidates_list

    when: params.star_fusion

    script:
    """
    STAR-Fusion \\
        --genome_lib_dir ${star_fusion_reference}\\
        --left_fq ${reads[0]} \\
        --right_fq ${reads[1]} \\
        --output_dir  $name
    """
}


/*
 * Fusion Catcher
*/
// Requires raw untrimmed files. FastQ files should not be merged!
process fusioncatcher {
    tag "$name"

    publishDir "${params.outdir}/FusionCatcher",  mode: 'copy'    

    input:
    set val (name), file(reads) from fusioncatcher_reads
    file fusioncatcher_data_dir from fusioncatcher_data_dir.collect()

    output:
    file '*.{txt,log,zip}' into fusioncatcher
    when: params.fusioncatcher

    script:
    if (params.singleEnd) {
    
    """
    mkdir ${reads}_data
    mv ${reads} ${reads}_data/
    fusioncatcher \\
        -d $fusioncatcher_data_dir \\
        -i ${reads}_data \\
        --threads ${task.cpus} \\
        -o $name \\
        --skip-blat \\
        --single-end \\
        ${params.fc_extra_options}
    """
    } else {
 
    """
    fusioncatcher \\
        -d $fusioncatcher_data_dir \\
        -i ${reads[0]},${reads[1]} \\
        --threads ${task.cpus} \\
        --${params.sensitivity} \\
        -o $name \\
        --skip-blat \\
        ${params.fc_extra_options}
    """
}
}



process fusion_genes_compare {
    
    publishDir "${params.outdir}/Comparative_shortlist",  mode: 'copy'
    
    input:
    file ('*star_fusion.fusion_candidates.final.abridged') from fusion_candidates_list.collect() 
    file ('*summary_candidate_fusions.txt') from fusioncatcher.collect()
    
    output:
    file '*fusion_comparison.txt'  
    
    when: params.star && params.inspector

    script:
    """
    fusion_genes_compare.py * 
    """

}


/*
 *  -  FusionInspector
 */
process fusioninspector {
    tag "$name"
    publishDir "${params.outdir}/FusionInspector",  mode: 'copy'  
    input:
    set val (name), file(reads) from fusion_inspector_reads
    file fusion_candidates
    file star_fusion_reference from star_fusion_reference_fusioninspector.collect()

    output:
    file '*' into fusioninspector_results

    when: params.inspector

    script:
    """
    FusionInspector \\
        --fusions $fusion_candidates \\
        --genome_lib $star_fusion_reference \\
        --left_fq ${reads[0]} \\
        --right_fq ${reads[1]} \\
        --out_dir . \\ 
        --out_prefix finspector \\
        --prep_for_IGV
    """
}


