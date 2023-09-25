process GTF_TO_REFFLAT {
    label 'process_low'

    conda "bioconda::ucsc-gtftogenepred=377"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-gtftogenepred:377--ha8a8165_5' :
        'quay.io/biocontainers/ucsc-gtftogenepred:377--ha8a8165_5' }"

    input:
    tuple val(meta), path (gtf)

    output:
    path('*.refflat'), emit: refflat

    script:
    def genepred = gtf + '.genepred'
    def refflat = gtf + '.refflat'
    """
    gtfToGenePred -genePredExt -geneNameAsName2 ${gtf} ${genepred}
    paste ${genepred} ${genepred} | cut -f12,16-25 > ${refflat}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtfToGenePred: 377
    END_VERSIONS
    """

    stub:
    def refflat = gtf + '.refflat'
    """
    touch ${refflat}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtfToGenePred: 377
    END_VERSIONS
    """
}
