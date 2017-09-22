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
params.outdir = './results'
params.email = false
params.star = true
params.inspector = false
params.fusioncatcher = true
params.sensitivity = 'sensitive'
params.clusterOptions = false
params.fc_extra_options = ''
Channel
    .fromFilePairs( params.reads, size: 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}" }
    .into { read_files_star_fusion; fusion_inspector_reads; fusioncatcher_reads}


/*
 * STAR-Fusion
 */
process star_fusion{
    tag "$name"
    publishDir 'Star_fusion',  mode: 'copy'
    
    input:
    set val (name), file(reads) from read_files_star_fusion

    output:
    file '*final.abridged*' into star_fusion_abridged
    file '*star-fusion.fusion_candidates.final.abridged.FFPM' into fusion_candidates

    when: params.star

    script:
    """
    STAR-Fusion \\
        --genome_lib_dir ${params.star_fusion_reference}\\
        --left_fq ${reads[0]} \\
        --right_fq ${reads[1]} \\
        --output_dir .
    for f in *
    do
    mv \$f $name\$f
    done
    """
}


/*
 *  -  FusionInspector
 */
process fusioninspector {
    tag "$name"
    publishDir 'FusionInspector',  mode: 'copy'  
    input:
    set val (name), file(reads) from fusion_inspector_reads
    file fusion_candidates

    output:
    file '*' into fusioninspector_results

    when: params.inspector

    script:
    """
    FusionInspector \\
        --fusions $fusion_candidates \\
        --genome_lib ${params.star_fusion_reference} \\
        --left_fq ${reads[0]} \\
        --right_fq ${reads[1]} \\
        --out_dir . \\ 
        --out_prefix finspector \\
        --prep_for_IGV
    """
}



/*
 * Fusion Catcher
*/
// Requires raw untrimmed files. FastQ files should not be merged!
process fusioncatcher {
    tag "$name"

    publishDir 'FuisonCatcher',  mode: 'copy'    

    input:
    set val (name), file(reads) from fusioncatcher_reads

    output:
    file '*.{txt,log,zip}' into fusioncatcher

    when: params.fusioncatcher

    script:
    """
    fusioncatcher \\
        -d ${params.fusioncatcher_data_dir} \\
        -i ${reads[0]},${reads[1]} \\
        --threads ${task.cpus} \\
        --${params.sensitivity} \\
        -o . \\
        ${params.fc_extra_options}
    for f in *.{txt,log,zip}
    do
    mv \$f $name\$f
    done
    """
}



