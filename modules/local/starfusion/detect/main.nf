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
    path reference

    output:
    tuple val(meta), path("*.fusion_predictions.tsv"), emit: fusions
    tuple val(meta), path("*.abridged.tsv")          , emit: abridged
    tuple val(meta), path("*.coding_effect.tsv")     , optional: true, emit: coding_effect
    path "*.version.txt"                             , emit: version

    script:
    def software    = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    STAR \\
	    --genomeDir $index \\
	    --readFilesIn $reads \\
	    --runThreadN $task.cpus \\
        $options.args

    STAR-Fusion \\
        --genome_lib_dir $reference \\
        $reads \\
        --CPU $task.cpus \\
        --output_dir . \\
        $options.args2

    mv star-fusion.fusion_predictions.tsv ${prefix}.starfusion.fusion_predictions.tsv
    mv star-fusion.fusion_predictions.abridged.tsv ${prefix}.starfusion.abridged.tsv
    mv star-fusion.fusion_predictions.abridged.coding_effect.tsv ${prefix}.starfusion.abridged.coding_effect.tsv
    
    echo \$(STAR-Fusion --version 2>&1) | grep -i 'version' | sed 's/STAR-Fusion version: //' > ${software}.version.txt
    """
}
