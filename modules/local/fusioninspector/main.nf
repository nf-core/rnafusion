process FUSIONINSPECTOR {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::dfam=3.3 bioconda::hmmer=3.3.2 bioconda::star-fusion=1.10.0 bioconda::trinity=date.2011_11_2 bioconda::samtools=1.9 bioconda::star=2.7.8a" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "docker.io/trinityctat/starfusion:1.10.1"
    } else {
        container "docker.io/trinityctat/starfusion:1.10.1"
    }

    input:
    tuple val(meta), path(reads)
    path fusion_list
    path index

    output:
    path "*"                , emit: output
    path "versions.yml"     , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta = params.single_end ? "--left_fq ${reads[0]}" : "--left_fq ${reads[0]} --right_fq ${reads[1]}"
    def args = task.ext.args ?: ''
    """
    FusionInspector \\
        --fusions ${fi_input_list} \\
        --genome_lib ${reference} \\
        $fasta \\
        --CPU ${task.cpus} \\
        -O . \\
        --out_prefix $prefix \\
        --vis $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR-Fusion: \$(echo STAR-Fusion --version 2>&1 | grep -i 'version' | sed 's/STAR-Fusion version: //')
    END_VERSIONS
    """
}
