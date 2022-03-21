process GTF_TO_REFFLAT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::ucsc-gtftogenepred=377" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-gtftogenepred:377--ha8a8165_5' :
        'quay.io/biocontainers/ucsc-gtftogenepred:377--ha8a8165_5' }"

    input:
    path gtf

    output:
    path('*.refflat'), emit: refflat

    script:
    def genepred = ${gtf}.getSimpleName() + '.genepred'
    def refflat = ${gtf}.getSimpleName() + '.refflat'
    """
    gtfToGenePred -genePredExt -geneNameAsName2 ${gtf} ${genepred}
    paste ${genepred} ${genepred} | cut -f12,16-25 > ${refflat}
    """
}
