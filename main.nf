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
params.star = false
params.fusioncatcher = true
Channel
    .fromFilePairs( params.reads, size: 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}" }
    .into { read_files_star_fusion; fusion_inspector_reads; fusioncatcher_reads}


/*
 * STEP 1 - STAR-Fusion
 */
process star_fusion{

    input:
    set val (name), file(read1), file(read2) from read_files_star_fusion

    output:
    file '*final.abridged*' into star_fusion_abridged
    file 'star-fusion.fusion_candidates.final.abridged.FFPM' into fusion_candidates

    when: params.star

    script:
    """
    STAR-Fusion \
        --genome_lib_dir $star_fusion_refrence \
        --left_fq $read1 \
        --right_fq $read2 \
        --output_dir $name
    """
}


/*
 *  -  FusionInspector
 */
process fusioninspector {

    input:
    set val (name), file(read1), file(read2) from fusion_inspector_reads
    file fusion_candidates

    output:
    file '*' into fusioninspector_results

    when: params.star

    script:
    """
    FusionInspector \
        --fusions $fusion_candidates \
        --genome_lib $star_fusion_refrence \
        --left_fq $read1 \
        --right_fq $read2 \
        --out_dir $my_FusionInspector_outdir \
        --out_prefix finspector \
        --prep_for_IGV
    """
}



/*
 * Fusion Catcher
*/
// Requires raw untrimmed files. FastQ files should not be merged!
process fusioncatcher {

    input:
    set val (name), file(read1), file(read2) from fusioncatcher_reads

    output:
    file '*.txt' into fusioncatcher

    when: params.fusioncatcher

    script:
    """
    fusioncatcher \
        -d $fusioncatcher_data_dir \
        -i ${read1},${read2} \
        --threads ${task.cpus} \
        --$SENSITIVITY \
        -o ${SAMPLE_NAME}/
    """
}



