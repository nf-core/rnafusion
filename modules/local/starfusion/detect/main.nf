process STARFUSION {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::dfam=3.7 bioconda::hmmer=3.4 bioconda::star-fusion=1.13.0 bioconda::trinity=2.15.1 bioconda::samtools=1.19.2 bioconda::star=2.7.11b"
    container 'docker.io/trinityctat/starfusion:1.13.0'

    input:
    tuple val(meta), path(reads), path(junction)
    path reference

    output:
    tuple val(meta), path("*.fusion_predictions.tsv")                   , emit: fusions
    tuple val(meta), path("*.abridged.tsv")                             , emit: abridged
    tuple val(meta), path("*.coding_effect.tsv")     , optional: true   , emit: coding_effect
    path "versions.yml"                                                 , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta = meta.single_end ? "--left_fq ${reads[0]}" : "--left_fq ${reads[0]} --right_fq ${reads[1]}"
    def args = task.ext.args ?: ''
    """
    STAR-Fusion \\
        --genome_lib_dir $reference \\
        $fasta \\
        -J $junction \\
        --CPU $task.cpus \\
        --examine_coding_effect \\
        --output_dir . \\
        $args

    mv star-fusion.fusion_predictions.tsv ${prefix}.starfusion.fusion_predictions.tsv
    mv star-fusion.fusion_predictions.abridged.tsv ${prefix}.starfusion.abridged.tsv
    mv star-fusion.fusion_predictions.abridged.coding_effect.tsv ${prefix}.starfusion.abridged.coding_effect.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR-Fusion: \$(STAR-Fusion --version 2>&1 | grep -i 'version' | sed 's/STAR-Fusion version: //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.starfusion.fusion_predictions.tsv
    touch ${prefix}.starfusion.abridged.tsv
    touch ${prefix}.starfusion.abridged.coding_effect.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR-Fusion: \$(STAR-Fusion --version 2>&1 | grep -i 'version' | sed 's/STAR-Fusion version: //')
    END_VERSIONS
    """
}


