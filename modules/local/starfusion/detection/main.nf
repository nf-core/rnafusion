// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process STARFUSION {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::star-fusion=1.10.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/star-fusion:1.10.0--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/star-fusion:1.10.0--hdfd78af_1"
    }

    input:
    tuple val(meta), path(reads)
    path genome_resource_lib

    output:
    tuple val(meta), path("*.fusion_predictions.tsv")   , emit: fusions
    path "*.version.txt"                                , emit: version

    script:
    def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def fastq       = meta.single_end ? "--left_fq ${reads[0]}" : "--left_fq ${reads[0]} --right_fq ${reads[1]}"

    """
    STAR-Fusion \\
        --genome_lib_dir $genome_resource_lib \\
        $fastq \\
        --output_dir . \\
        $options.args

    echo \$(STAR-Fusion --version 2>&1) | grep -i 'version' | sed 's/STAR-Fusion version: //' > ${software}.version.txt
    """
}
