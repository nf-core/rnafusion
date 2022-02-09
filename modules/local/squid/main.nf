
process SQUID {
    tag "squid"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::squid=1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/squid:1.5--haac9b31_4' :
        'https://depot.galaxyproject.org/singularity/squid:1.5--haac9b31_4' }"



    input:
    tuple val(meta), path(bam)
    path gtf

    output:
    tuple val(meta), path("*_sv.txt"), emit: fusions
    tuple val(meta), path("*_annotated.txt"), emit: fusions_annotated
    path  "versions.yml"          , emit: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    squid -b $bam -o ${prefix}_fusions
    AnnotateSQUIDOutput.py $gtf ${prefix}_fusions_sv.txt ${prefix}_fusions_annotated.txt
    """
}
