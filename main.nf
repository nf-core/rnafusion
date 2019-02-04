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

nfcore_logo = """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

nf-core/rnafusion v${workflow.manifest.version}
======================================================="""

def helpMessage() {
    nfcore_help = """
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/rnafusion --reads '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: standard, conda, docker, singularity, awsbatch, test

    Tool flags:
      --star_fusion                 Run STAR-Fusion
      --fusioncatcher               Run FusionCatcher
        --fc_extra_options          Extra parameters for FusionCatcher. Can be found at https://github.com/ndaniel/fusioncatcher/blob/master/doc/manual.md
      --ericscript                  Run Ericscript
      --pizzly                      Run Pizzly
      --squid                       Run Squid
      --test                        Runs only specific fusion tool/s and not the whole pipeline. Only works on tool flags.
      --tool_cutoff                 Number of tools to pass threshold required for showing fusion in the report. [Default = 2]

    Visualization flags:
      --fusion_inspector            Run Fusion-Inspector

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta                       Path to Fasta reference
      --gtf                         Path to GTF annotation
      --star_index                  Path to STAR-Index reference
      --star_fusion_ref             Path to STAR-Fusion reference
      --fusioncatcher_ref           Path to Fusioncatcher reference
      --ericscript_ref              Path to Ericscript reference
      --pizzly_fasta                Path to Pizzly FASTA reference
      --pizzly_gtf                  Path to Pizzly GTF annotation

    Options:
      --genome                      Name of iGenomes reference
      --read_length                 Length of the reads. Default: 100
      --singleEnd                   Specifies that the input is single end reads

    Other Options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
    log.info "${nfcore_logo}${nfcore_help}"
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
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.running_tools = []
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false

multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")

// Reference variables required by tools
// These are needed in order to run the pipeline
pizzly_fasta = false
pizzly_gtf = false
star_fusion_ref = false
fusioncatcher_ref = false
fusion_inspector_ref = false
ericscript_ref = false

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
if (params.fasta) {
    fasta = file(params.fasta)
    if(!fasta.exists()) exit 1, "Fasta file not found: ${params.fasta}"
} else {
    if (!params.genome) exit 1, "You have to specify either fasta file or Genome version!"
}

if (params.gtf) {
    gtf = file(params.gtf)
    if(!gtf.exists()) exit 1, "GTF file not found: ${params.fasta}"
}

if (params.star_fusion) {
    params.running_tools.add("STAR-Fusion")
    if (!params.star_fusion_ref) {
        exit 1, "Star-Fusion reference not specified!"
    } else {
        star_fusion_ref = Channel
            .fromPath(params.star_fusion_ref)
            .ifEmpty { exit 1, "Star-Fusion reference directory not found!" }
    }
}

if (params.fusioncatcher) {
    params.running_tools.add("Fusioncatcher")
    if (!params.fusioncatcher_ref) {
        exit 1, "Fusioncatcher data directory not specified!"
    } else {
        fusioncatcher_ref = Channel
            .fromPath(params.fusioncatcher_ref)
            .ifEmpty { exit 1, "Fusioncatcher data directory not found!" }
    }
}

if (params.ericscript) {
    params.running_tools.add("Ericscript")
    if (!params.ericscript_ref) {
        exit 1, "Reference not specified!"
    } else {
        ericscript_ref = Channel
            .fromPath(params.ericscript_ref)
            .ifEmpty { exit 1, "Ericscript reference not found" }
    }
}

if (params.pizzly) {
    params.running_tools.add("Pizzly")
    if (params.pizzly_fasta) {
        pizzly_fasta = Channel
            .fromPath(params.pizzly_fasta)
            .ifEmpty { exit 1, "Pizzly FASTA file not found!" }
    }

    if (params.pizzly_gtf) {
        pizzly_gtf = Channel
            .fromPath(params.pizzly_gtf)
            .ifEmpty { exit 1, "Pizzly GTF file not found!" }
    }
}

if (params.squid) {
    params.running_tools.add("Squid")
}

if (params.fusion_inspector) {
    params.running_tools.add("FusionInspector")
    if (!params.star_fusion_ref) {
        exit 1, "Reference not specified (using star-fusion reference path)!"
    } else {
        fusion_inspector_ref = Channel
            .fromPath(params.star_fusion_ref)
            .ifEmpty { exit 1, "Fusion-Inspector reference not found" }
    }
}

/*
 * Create a channel for input read files
 */
Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { read_files_fastqc; read_files_summary; read_files_multiqc; read_files_star_fusion; read_files_fusioncatcher; 
            read_files_gfusion; read_files_fusion_inspector; read_files_ericscript; read_files_pizzly; read_files_squid }


// Header log info
log.info nfcore_logo
def summary = [:]
summary['Pipeline Name']  = 'nf-core/rnafusion'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Fasta Ref']    = params.fasta
summary['GTF Ref']      = params.gtf
summary['STAR Index']   = params.star_index ? params.star_index : 'Not specified, building'
summary['Tools']        = params.running_tools.size() == 0 ? 'None' : params.running_tools.join(", ")
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
 * PREPROCESSING
 ************************************************************/

/*
 * Build STAR index
 */
if (params.star_index) {
    Channel
        .fromPath(params.star_index)
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
        .into { star_index_squid; star_index_star_fusion }
} else {
    process build_star_index {
        tag "$fasta"
        publishDir "${params.outdir}/star_index", mode: 'copy'

        input:
        file fasta
        file gtf

        output:
        file "star" into star_index_squid, star_index_star_fusion

        script:
        def avail_mem = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
        """
        mkdir star
        STAR \\
            --runMode genomeGenerate \\
            --runThreadN ${task.cpus} \\
            --sjdbGTFfile ${gtf} \\
            --sjdbOverhang ${params.read_length - 1} \\
            --genomeDir star/ \\
            --genomeFastaFiles ${fasta} \\
            $avail_mem
        """
    }
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
    file star_index_star_fusion
    file reference from star_fusion_ref

    output:
    file '*fusion_predictions.tsv' into star_fusion_fusions
    file '*.{tsv,txt}' into star_fusion_output

    script:
    def avail_mem = task.memory ? "--limitBAMsortRAM ${task.memory.toBytes() - 100000000}" : ''
    option = params.singleEnd ? "--left_fq ${reads[0]}" : "--left_fq ${reads[0]} --right_fq ${reads[1]}"
    """
    STAR \\
        --genomeDir ${star_index_star_fusion} \\
        --readFilesIn ${reads} \\
        --twopassMode Basic \\
        --outReadsUnmapped None \\
        --chimSegmentMin 12 \\
        --chimJunctionOverhangMin 12 \\
        --alignSJDBoverhangMin 10 \\
        --alignMatesGapMax 100000 \\
        --alignIntronMax 100000 \\
        --chimSegmentReadGapMax 3 \\
        --alignSJstitchMismatchNmax 5 -1 5 5 \\
        --runThreadN ${task.cpus} \\
        --outSAMstrandField intronMotif ${avail_mem} \\
        --readFilesCommand zcat \\
        --chimOutJunctionFormat 1

    STAR-Fusion \\
        --genome_lib_dir ${reference} \\
        -J Chimeric.out.junction \\
        ${option} \\
        --CPU ${task.cpus} \\
        --examine_coding_effect \\
        --output_dir .
    """
}

/*
 * Fusioncatcher
 */
process fusioncatcher {
    tag "$name"
    publishDir "${params.outdir}/tools/Fusioncatcher", mode: 'copy'

    when:
    params.fusioncatcher || (params.fusioncatcher && params.test)

    input:
    set val(name), file(reads) from read_files_fusioncatcher
    file data_dir from fusioncatcher_ref

    output:
    file 'final-list_candidate-fusion-genes.txt' into fusioncatcher_fusions
    file '*.{txt,zip,log}' into fusioncatcher_output

    script:
    option = params.singleEnd ? reads[0] : "${reads[0]},${reads[1]}"
    """
    fusioncatcher \\
        -d ${data_dir} \\
        -i ${option} \\
        --threads ${task.cpus} \\
        -o . \\
        --skip-blat \\
        ${params.fc_extra_options}
    """
}

/*
 * Ericscript
 */
process ericscript {
    tag "$name"
    publishDir "${params.outdir}/tools/Ericscript", mode: 'copy'

    when:
    params.ericscript && (!params.singleEnd || params.test)

    input:
    set val(name), file(reads) from read_files_ericscript
    file reference from ericscript_ref

    output:
    file './tmp/fusions.results.filtered.tsv' into ericscript_fusions
    file './tmp/fusions.results.total.tsv' into ericscript_output

    script:
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

/*
 * Pizzly
 */
process pizzly {
    tag "$name"
    publishDir "${params.outdir}/tools/Pizzly", mode: 'copy'

    when:
    params.pizzly && (!params.singleEnd || params.test)

    input:
    set val(name), file(reads) from read_files_pizzly
    file fasta from pizzly_fasta
    file gtf from pizzly_gtf
    
    output:
    file 'pizzly_fusions.json' into pizzly_fusions
    file '*.{json,txt}' into pizzly_output

    script:
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

/*
 * Squid
 */
process squid {
    tag "$name"
    publishDir "${params.outdir}/tools/Squid", mode: 'copy'

    when:
    params.squid && (!params.singleEnd || params.test)

    input:
    set val(name), file(reads) from read_files_squid
    file star_index_squid
    file gtf
    
    output:
    file '*_annotated.txt' into squid_fusions
    file '*.txt' into squid_output

    script:
    def avail_mem = task.memory ? "--limitBAMsortRAM ${task.memory.toBytes() - 100000000}" : ''
    """
    STAR \\
        --genomeDir ${star_index_squid} \\
        --sjdbGTFfile ${gtf} \\
        --runThreadN ${task.cpus} \\
        --readFilesIn ${reads[0]} ${reads[1]} \\
        --twopassMode Basic \\
        --chimOutType SeparateSAMold --chimSegmentMin 20 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 --outReadsUnmapped Fastx --outSAMstrandField intronMotif \\
        --outSAMtype BAM SortedByCoordinate ${avail_mem} \\
        --readFilesCommand zcat
    mv Aligned.sortedByCoord.out.bam ${name}Aligned.sortedByCoord.out.bam
    samtools view -bS Chimeric.out.sam > ${name}Chimeric.out.bam
    squid -b ${name}Aligned.sortedByCoord.out.bam -c ${name}Chimeric.out.bam -o fusions
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
    extra = params.tool_cutoff ? "-t ${params.tool_cutoff}" : "" 
    """
    transformer.py -i ${fusioncatcher} -t fusioncatcher
    transformer.py -i ${star_fusion} -t star_fusion
    transformer.py -i ${ericscript} -t ericscript
    transformer.py -i ${pizzly} -t pizzly
    transformer.py -i ${squid} -t squid
    generate_report.py fusions.txt summary.yaml -s ${name} -o . ${extra}
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
    params.fusion_inspector && (!params.singleEnd || params.test)

    input:
    set val(name), file(reads) from read_files_fusion_inspector
    file reference from fusion_inspector_ref
    file summary_fusions

    output:
    file '*.{fa,gtf,bed,bam,bai,txt}' into fusion_inspector_output

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