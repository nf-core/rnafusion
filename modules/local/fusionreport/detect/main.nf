process FUSIONREPORT {
    tag 'fusionreport'
    label 'process_medium'

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda (params.enable_conda ? 'bioconda::star=2.7.9a' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/rannickscilifelab/fusion-report:2.1.5' :
        'docker.io/rannickscilifelab/fusion-report:2.1.5' }"


    input:
    tuple val(meta), path(reads)
    path(fusionreport_ref)
    file(arriba_tools)

    output:
    path "*"                , emit: reference
    path "versions.yml"     , emit: versions

    script:

    """
    fusion_report run $meta.id . $fusionreport_ref --arriba $arriba_tools --allow-multiple-gene-symbols

    fusion_report --version > versions.yml
    """
}
