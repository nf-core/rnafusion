process MULTIQC {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publishDir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publishDir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? 'bioconda::multiqc=1.11' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/multiqc:1.11--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0"
    }

    input:
    path multiqc_files

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    script:
    """
    multiqc -f $options.args .
    multiqc --version | sed -e "s/multiqc, version //g" > versions.yml
    """
}
