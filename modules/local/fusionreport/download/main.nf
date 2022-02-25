process FUSIONREPORT_DOWNLOAD {
    tag 'fusionreport'
    label 'process_medium'

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda (params.enable_conda ? 'bioconda::star=2.7.9a' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/rannickscilifelab/fusion-report:2.1.5' :
        'docker.io/rannickscilifelab/fusion-report:2.1.5' }"


    input:
    val(username)
    val(passwd)

    output:
    path "*"                , emit: reference
    path "versions.yml"     , emit: versions

    script:
    """
    fusion_report download --cosmic_usr $username --cosmic_passwd $passwd .

    fusion_report --version > versions.yml
    """
}
