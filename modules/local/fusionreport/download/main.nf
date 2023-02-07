process FUSIONREPORT_DOWNLOAD {
    tag 'fusionreport'
    label 'process_medium'

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda "bioconda::star=2.7.9a"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/clinicalgenomics/fusion-report:2.1.5patched' :
        'docker.io/clinicalgenomics/fusion-report:2.1.5patched' }"


    input:
    val(username)
    val(passwd)

    output:
    path "*"                , emit: reference
    path "versions.yml"     , emit: versions

    script:
    """
    fusion_report download --cosmic_usr $username --cosmic_passwd $passwd .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusion_report: \$(fusion_report --version | sed 's/fusion-report //')
    END_VERSIONS
    """
}
