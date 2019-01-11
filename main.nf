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
      --singleEnd                   Specifies that the input is single end reads
      --star_fusion                 [bool] Run STAR-Fusion. Default: False
      --fusioncatcher               [bool] Run FusionCatcher. Default: False
        --fc_extra_options          Extra parameters for FusionCatcher. Can be found at https://github.com/ndaniel/fusioncatcher/blob/master/doc/manual.md
      --fusion_inspector            [bool] Run Fusion-Inspectro. Default: False
      --ericscript                  [bool] Run Ericscript. Default: False
      --pizzly                      [Bool] Run Pizzly. Default: False
      --squid                       [Bool] Run Squid. Default: False
      --test                        [bool] Run in test mode

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta                       Path to Fasta reference

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

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
params.test = false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false

multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")

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

if (params.fasta) {
    fasta = file(params.fasta)
    if(!fasta.exists()) exit 1, "Fasta file not found: ${params.fasta}"
}

if (params.gtf) {
    gtf = file(params.gtf)
    if(!gtf.exists()) exit 1, "GTF file not found: ${params.fasta}"
}

if(params.star_index){
    star_index = Channel
        .fromPath(params.star_index)
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
}

star_fusion_ref = ''
if (params.star_fusion) {
    if (!params.star_fusion_ref) {
        exit 1, "Star-Fusion reference not specified!"
    } else {
        star_fusion_ref = Channel
            .fromPath(params.star_fusion_ref)
            .ifEmpty { exit 1, "Stat-Fusion reference directory not found!" }
    }
}

fusioncatcher_ref = ''
if (params.fusioncatcher) {
    if (!params.fusioncatcher_ref) {
        exit 1, "Fusion catcher data directory not specified!"
    } else {
        fusioncatcher_ref = Channel
            .fromPath(params.fusioncatcher_ref)
            .ifEmpty { exit 1, "Fusion catcher data directory not found!" }
    }
}

fusion_inspector_ref = ''
if (params.fusion_inspector) {
    if (!params.star_fusion_ref) {
        exit 1, "Reference not specified (using star-fusion reference path)!"
    } else {
        fusion_inspector_ref = Channel
            .fromPath(params.star_fusion_ref)
            .ifEmpty { exit 1, "Fusion-Inspector reference not found" }
    }
}

ericscript_ref = ''
if (params.ericscript) {
    if (!params.ericscript_ref) {
        exit 1, "Reference not specified!"
    } else {
        ericscript_ref = Channel
            .fromPath(params.ericscript_ref)
            .ifEmpty { exit 1, "Ericscript reference not found" }
    }
}

pizzly_fasta = ''
pizzly_gtf = ''
if (params.pizzly) {
    if (!params.pizzly_ref) {
        exit 1, "No pizzly reference folder defined"
    } else {
        pizzly_fasta = Channel
            .fromPath("$params.pizzly_ref/Homo_sapiens.GRCh38.cdna.all.fa.gz")
            .ifEmpty { exit 1, "FASTA file not found!" }
        pizzly_gtf = Channel
            .fromPath("$params.pizzly_ref/Homo_sapiens.GRCh38.94.gtf")
            .ifEmpty { exit 1, "GTF file not found!" }
    }
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
             .into { read_files_fastqc; read_files_summary; read_files_multiqc; read_files_star_fusion; read_files_fusioncatcher; 
                     read_files_gfusion; read_files_fusion_inspector; read_files_ericscript; read_files_pizzly; read_files_squid }
     } else {
         Channel
             .from(params.readPaths)
             .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
             .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
             .into { read_files_fastqc; read_files_summary; read_files_multiqc; read_files_star_fusion; read_files_fusioncatcher; 
                     read_files_gfusion; read_files_fusion_inspector; read_files_ericscript; read_files_pizzly; read_files_squid }
     }
 } else {
     Channel
         .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
         .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
         .into { read_files_fastqc; read_files_summary; read_files_multiqc; read_files_star_fusion; read_files_fusioncatcher; 
                     read_files_gfusion; read_files_fusion_inspector; read_files_ericscript; read_files_pizzly; read_files_squid }
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
    publishDir "${params.outdir}/tools/StarFusion", mode: 'copy'

    when:
    params.star_fusion || (params.star_fusion && params.test)

    input:
    set val(name), file(reads) from read_files_star_fusion
    file reference from star_fusion_ref.collect()

    output:
    file '*fusion_predictions.tsv' into star_fusion_fusions
    file '*' into star_fusion_output

    script:
    if (params.singleEnd) {
        """
        STAR-Fusion \\
            --genome_lib_dir ${reference} \\
            --left_fq ${reads[0]} \\
            --CPU  ${task.cpus} \\
            --examine_coding_effect \\
            --output_dir .
        """
    } else {
        """
        STAR-Fusion \\
            --genome_lib_dir ${reference} \\
            --left_fq ${reads[0]} \\
            --right_fq ${reads[1]} \\
            --CPU  ${task.cpus} \\
            --examine_coding_effect \\
            --output_dir .
        """
    } 
}

/*
 * Fusion Catcher
 */
process fusioncatcher {
    tag "$name"
    publishDir "${params.outdir}/tools/Fusioncatcher", mode: 'copy'

    when:
    params.fusioncatcher || (params.fusioncatcher && params.test)

    input:
    set val(name), file(reads) from read_files_fusioncatcher
    file data_dir from fusioncatcher_ref.collect()

    output:
    file 'final-list_candidate-fusion-genes.txt' into fusioncatcher_fusions
    file '*' into fusioncatcher_output

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
 * Ericscript
 */
process ericscript {
    tag "$name"
    publishDir "${params.outdir}/tools/Ericscript", mode: 'copy'

    when:
    params.ericscript || (params.ericscript && params.test)

    input:
    set val(name), file(reads) from read_files_ericscript
    file reference from ericscript_ref.collect()

    output:
    file './tmp/fusions.results.total.tsv' into ericscript_fusions
    file '*' into ericscript_output

    script:
    if (!params.singleEnd) {
        """
        ericscript.pl \\
            -db ${reference} \\
            -name fusions \\
            -p ${task.cpus} \\
            -o ./tmp \\
            ${reads[0]} \\
            ${reads[1]}
        """
    }
}

/*
 * Pizzly
 */
process pizzly {
    tag "$name"
    publishDir "${params.outdir}/tools/Pizzly", mode: 'copy'

    when:
    params.pizzly || (params.pizzly && params.test)

    input:
    set val(name), file(reads) from read_files_pizzly
    file fasta from pizzly_fasta.collect()
    file gtf from pizzly_gtf.collect()
    
    output:
    file '*.unfiltered.json' into pizzly_fusions
    file '*.{json,txt,tsv,fasta}' into pizzly_output

    script:
    if (!params.singleEnd) {
        """
        kallisto index -i index.idx -k ${params.pizzly_k} ${fasta}
        kallisto quant -t ${task.cpus} -i index.idx --fusion -o output ${reads[0]} ${reads[1]}
        pizzly -k ${params.pizzly_k} \\
            --gtf ${gtf} \\
            --cache index.cache.txt \\
            --align-score 2 \\
            --insert-size 400 \\
            --fasta ${fasta} \\
            --output pizzly_fusions output/fusion.txt
        """
    }
}

/*
 * Squid
 */
process squid {
    tag "$name"
    publishDir "${params.outdir}/tools/Squid", mode: 'copy'

    when:
    params.squid || (params.squid && params.test)

    input:
    set val(name), file(reads) from read_files_squid
    file index from star_index.collect()
    file gtf
    
    output:
    file '*_annotated.txt' into squid_fusions
    file '*' into squid_output

    script:
    def avail_mem = task.memory ? "--limitBAMsortRAM ${task.memory.toBytes() - 100000000}" : ''
    """
    STAR \\
        --genomeDir ${index} \\
        --sjdbGTFfile ${gtf} \\
        --runThreadN ${task.cpus} \\
        --readFilesIn ${reads[0]} ${reads[1]} \\
        --twopassMode Basic \\
        --chimOutType SeparateSAMold --chimSegmentMin 20 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 --outReadsUnmapped Fastx --outSAMstrandField intronMotif \\
        --outSAMtype BAM Unsorted ${avail_mem} \\
        --readFilesCommand zcat
    samtools sort Aligned.out.bam > Aligned.out.sorted.bam
    samtools view -Shb Chimeric.out.sam > Chimeric.out.bam
    squid -b Aligned.out.sorted.bam -c Chimeric.out.bam -o fusions
    AnnotateSQUIDOutput.py ${gtf} fusions_sv.txt fusions_annotated.txt
    """
}

/*************************************************************
 * Summarizing results from tools
 ************************************************************/
process summary {
    tag "$name"
    publishDir "${params.outdir}/Report-${name}", mode: 'copy'
 
    when:
    !params.test && (params.fusioncatcher || params.star_fusion || params.ericscript || params.pizzly || params.squid)
    
    input:
    set val(name), file(reads) from read_files_summary
    file fusioncatcher from fusioncatcher_fusions.ifEmpty('')
    file star_fusion from star_fusion_fusions.ifEmpty('')
    file ericscript from ericscript_fusions.ifEmpty('')
    file pizzly from pizzly_fusions.ifEmpty('')
    file squid from squid_fusions.ifEmpty('')

    output:
    file 'fusions.txt' into summary_fusions
    file 'summary.yaml' into summary_fusions_mq
    file '*.html' into report
    
    script:
    """
    transformer.py -i ${fusioncatcher} -t fusioncatcher
    transformer.py -i ${star_fusion} -t star_fusion
    transformer.py -i ${ericscript} -t ericscript
    transformer.py -i ${pizzly} -t pizzly
    transformer.py -i ${squid} -t squid
    generate_report.py fusions.txt summary.yaml -s ${name} -o .
    """
}

/*************************************************************
 * Visualization
 ************************************************************/

/*
 * Fusion Inspector
 */
process fusion_inspector {
    tag "$name"
    publishDir "${params.outdir}/tools/FusionInspector", mode: 'copy'

    when:
    params.fusion_inspector || (params.fusion_inspector && params.test)

    input:
    set val(name), file(reads) from read_files_fusion_inspector
    file reference from fusion_inspector_ref.collect()
    file summary_fusions

    output:
    file '*' into fusion_inspector_output

    script:
    """
    FusionInspector \\
        --fusions ${summary_fusions} \\
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

    when:
    !params.test

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    cat $baseDir/tools/fusioncatcher/Dockerfile | grep "VERSION" > v_fusioncatcher.txt
    cat $baseDir/tools/fusion-inspector/environment.yml | grep "fusion-inspector" > v_fusion_inspector.txt
    cat $baseDir/tools/star-fusion/environment.yml | grep "star-fusion" > v_star_fusion.txt
    cat $baseDir/tools/ericscript/environment.yml | grep "ericscript" > v_ericscript.txt
    cat $baseDir/tools/pizzly/environment.yml | grep "pizzly" > v_pizzly.txt
    cat $baseDir/tools/squid/environment.yml | grep "squid" > v_squid.txt
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
    !params.test

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
 * MultiQC
 */
process multiqc {
    tag "$name"
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    when:
    !params.test

    input:
    set val(name), file(reads) from read_files_multiqc
    file multiqc_config
    file ('fastqc/*') from fastqc_results.collect()
    file ('software_versions/*') from software_versions_yaml
    file workflow_summary from create_workflow_summary(summary)
    file fusions_mq from summary_fusions_mq.ifEmpty('')

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    create_mqc_section.py -i ${fusions_mq} -s "${name}"
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}

/*
 * Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    when:
    !params.test

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