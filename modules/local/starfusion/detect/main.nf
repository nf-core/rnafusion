process STARFUSION {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publishDir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publishDir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::dfam=3.3 bioconda::hmmer=3.3.2 bioconda::star-fusion=1.10.0 bioconda::trinity=date.2011_11_2 bioconda::samtools=1.9 bioconda::star=2.7.8a" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-57582e8bdbf51679bdcff9de91ae016a44e322de:b978a3a1e8715581329a267a2c9904574384180a-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-57582e8bdbf51679bdcff9de91ae016a44e322de:b978a3a1e8715581329a267a2c9904574384180a-0"
    }

    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path("*.fusion_predictions.tsv"), emit: fusions
    tuple val(meta), path("*.abridged.tsv")          , emit: abridged
    tuple val(meta), path("*.coding_effect.tsv")     , optional: true, emit: coding_effect
    path "versions.yml"                              , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    STAR \\
	    --genomeDir $index \\
	    --readFilesIn $reads \\
	    --runThreadN $task.cpus \\
        $args

    STAR-Fusion \\
        --genome_lib_dir $reference \\
        $reads \\
        --CPU $task.cpus \\
        --output_dir . \\
        $args2

    mv star-fusion.fusion_predictions.tsv ${prefix}.starfusion.fusion_predictions.tsv
    mv star-fusion.fusion_predictions.abridged.tsv ${prefix}.starfusion.abridged.tsv
    mv star-fusion.fusion_predictions.abridged.coding_effect.tsv ${prefix}.starfusion.abridged.coding_effect.tsv

    echo \$(STAR-Fusion --version 2>&1) | grep -i 'version' | sed 's/STAR-Fusion version: //' > versions.yml
    """
}
