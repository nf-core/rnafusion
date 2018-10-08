#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/rnafusion
========================================================================================
 nf-core/rnafusion Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/rnafusion
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     nf-core/rnafusion v${workflow.manifest.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/rnafusion --reads '*_R{1,2}.fastq.gz' -profile standard,docker

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --genome                      Name of iGenomes reference
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: standard, conda, docker, singularity, awsbatch, test

    Options:
      --all                         [bool] Run all tools
      --singleEnd                   Specifies that the input is single end reads
      --star_fusion                 [bool] Run STAR-Fusion. Default: False
      --fusioncatcher               [bool] Run FusionCatcher. Default: False
        --fc_extra_options          Extra parameters for FusionCatcher. Can be found at https://github.com/ndaniel/fusioncatcher/blob/master/doc/manual.md
      --fusion_inspector            [bool] Run Fusion-Inspectro. Default: False

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta                       Path to Fasta reference

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    QC options:
      --skip_qc                     Skip all QC including MultiQC

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false
params.all = false
params.skip_qc = false
params.star_fusion = false
params.fusioncatcher = false
params.fusion_inspector = false
params.fc_extra_options = ''

multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")

// Validate inputs
if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}
// AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Check workDir/outdir paths to be S3 buckets if running on AWSBatch
// related: https://github.com/nextflow-io/nextflow/issues/813
if( workflow.profile == 'awsbatch') {
    if(!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}

// Validate pipeline variables
// These variable have to be defined in the profile configuration which is referenced in nextflow.config
if (!params.genome) {
    exit 1, "Mandatory parameter genome not specified!"    
}

if (!params.star_fusion_ref) {
    exit 1, "Star-Fusion reference not specified!"
} else {
     Channel
        .fromPath(params.star_fusion_ref)
        .ifEmpty { exit 1, "Stat-Fusion reference directory not found!" }
        .into { star_fusion_ref; fusion_inspector_ref }
}

if (!params.fusioncatcher_dir) {
    exit 1, "Fusion catcher data directory not specified!"
} else {
    fusioncatcher_dir = Channel
        .fromPath(params.fusioncatcher_dir)
        .ifEmpty { exit 1, "Fusion catcher data directory not found!" }
}

/*
 * Create a channel for input read files
 */
 if(params.readPaths){
     if(params.singleEnd){
         Channel
             .from(params.readPaths)
             .map { row -> [ row[0], [file(row[1][0])]] }
             .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
             .into { read_files_fastqc; read_files_star_fusion; read_files_fusioncatcher; read_files_gfusion; read_files_fusion_inspector }
     } else {
         Channel
             .from(params.readPaths)
             .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
             .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
             .into { read_files_fastqc; read_files_star_fusion; read_files_fusioncatcher; read_files_gfusion; read_files_fusion_inspector }
     }
 } else {
     Channel
         .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
         .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
         .into { read_files_fastqc; read_files_star_fusion; read_files_fusioncatcher; read_files_gfusion; read_files_fusion_inspector }
 }


// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

nf-core/rnafusion v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'nf-core/rnafusion'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Fasta Ref']    = params.fasta
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


def create_workflow_summary(summary) {

    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-rnafusion-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/rnafusion Workflow Summary'
    section_href: 'https://github.com/nf-core/rnafusion'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}

/*************************************************************
 * Fusion pipeline
 ************************************************************/

/*
 * STAR-Fusion
 */
process star_fusion {
    tag "$name"
    publishDir "${params.outdir}/Star_fusion", mode: 'copy'

    when:
    params.star_fusion || params.all

    input:
    set val(name), file(reads) from read_files_star_fusion
    file reference from star_fusion_ref.collect()

    output:
    file '*fusion_predictions.tsv' into star_fusion_predictions
    file '*fusion_predictions.abridged.tsv' into star_fusion_abridged

    script:
    if (params.singleEnd) {
        """
        STAR-Fusion \\
        --genome_lib_dir ${reference} \\
        --left_fq ${reads[0]} \\
        --CPU  ${task.cpus} \\
        --output_dir .
        """
    } else {
        """
        STAR-Fusion \\
        --genome_lib_dir ${reference} \\
        --left_fq ${reads[0]} \\
        --right_fq ${reads[1]} \\
        --CPU  ${task.cpus} \\
        --output_dir .
        """
    } 
}

/*
 * Fusion Catcher
 */
process fusioncatcher {
    tag "$name"
    publishDir "${params.outdir}/Fusioncatcher", mode: 'copy'

    when:
    params.fusioncatcher || params.all

    input:
    set val(name), file(reads) from read_files_fusioncatcher
    file data_dir from fusioncatcher_dir.collect()

    output:
    file 'final-list_candidate-fusion-genes.txt' into fusioncatcher_candidates
    file '*.{txt,log,zip}' into fusioncatcher_output

    script:
    if (params.singleEnd) {
        """
        fusioncatcher \\
        -d ${data_dir} \\
        -i ${reads[0]} \\
        --threads ${task.cpus} \\
        -o . \\
        --skip-blat \\
        --single-end \\
        ${params.fc_extra_options}
        """
    } else {
        """
        fusioncatcher \\
        -d ${data_dir} \\
        -i ${reads[0]},${reads[1]} \\
        --threads ${task.cpus} \\
        -o . \\
        --skip-blat \\
        ${params.fc_extra_options}
        """
    }
}

/*
 * Fusion Inspector preprocess
 */
process fusion_inspector_preprocess {
    tag "$name"
    publishDir "${params.outdir}/Transformers", mode: 'copy'

    when:
    params.fusion_inspector || params.all

    input:
    file fusioncatcher_candidates
    file star_fusion_abridged

    output:
    file 'fusions.txt' into fusions
    file 'summary.yaml' into fusions_summary
    
    script:
    """
    transformer.py -i ${fusioncatcher_candidates} -t fusioncatcher
    transformer.py -i ${star_fusion_abridged} -t star_fusion
    """
}

/*
 * Fusion Inspector
 */
process fusion_inspector {
    tag "$name"
    publishDir "${params.outdir}/FusionInspector", mode: 'copy'

    when:
    params.fusion_inspector || params.all

    input:
    set val(name), file(reads) from read_files_fusion_inspector
    file reference from fusion_inspector_ref
    file fusions

    output:
    file '*' into fusion_inspector_output

    script:
    """
    FusionInspector \\
        --fusions ${fusions} \\
        --genome_lib ${reference} \\
        --left_fq ${reads[0]} \\
        --right_fq ${reads[1]} \\
        --CPU ${task.cpus} \\
        --out_dir . \\
        --out_prefix finspector \\
        --prep_for_IGV
    """
}

/*************************************************************
 * Building report
 ************************************************************/

/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    STAR-Fusion --version > v_star_fusion.txt
    fusioncatcher --version > v_fusioncatcher.txt
    cat /environment.yml | grep 'fusion-inspector' > v_fusion_inspector.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}

/*
 * FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    when:
    !params.skip_qc || params.all

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}

/*
 * Fusion gene compare
 * Builds MultiQC custom section
 */
process fusion_gene_compare {
    publishDir "${params.outdir}/FusionGeneCompare", mode: 'copy'

    when:
    params.fusion_inspector || params.all

    input:
    file fusions_summary

    output:
    file 'fusion_genes_mqc.yaml' into fusion_genes_mqc

    script:
    """
    fusion_genes_compare.py -i ${fusions_summary} -s test_sample
    """
}

/*
 * MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    when:
    !params.skip_qc || params.all

    input:
    file multiqc_config
    file ('fastqc/*') from fastqc_results.collect()
    file ('software_versions/*') from software_versions_yaml
    file workflow_summary from create_workflow_summary(summary)
    file fusion_genes_mqc

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}

/*
 * Output Description HTML
 */
process output_documentation {
    tag "$prefix"
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    file output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/rnafusion] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/rnafusion] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/rnafusion] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/rnafusion] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[nf-core/rnafusion] Pipeline Complete"

}