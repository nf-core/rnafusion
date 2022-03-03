process FUSIONREPORT {
    tag "$meta.id"
    label 'process_medium'

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda (params.enable_conda ? 'bioconda::star=2.7.9a' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/rannickscilifelab/fusion-report:2.1.5' :
        'docker.io/rannickscilifelab/fusion-report:2.1.5' }"


    input:
    tuple val(meta), path(reads)
    path(fusionreport_ref)
    path(arriba_fusions)
    path(pizzly_fusions)
    path(squid_fusions)
    path(starfusion_fusions)
    path(fusioncatcher_fusions)

    output:
    path "versions.yml"                                  , emit: versions
    tuple val(meta), path("*fusionreport.tsv")           , emit: fusion_list
    tuple val(meta), path("*fusionreport_filtered.tsv")  , emit: fusion_list_filtered
    tuple val(meta), path("*.html")                      , emit: report

    script:
    def tools = params.arriba         ? "--arriba ${arriba_fusions} " : ''
    tools    += params.pizzly         ? "--pizzly ${pizzly_fusions} " : ''
    tools    += params.squid          ? "--squid ${squid_fusions} " : ''
    tools    += params.starfusion     ? "--starfusion ${starfusion_fusions} " : ''
    tools    += params.fusioncatcher  ? "--fusioncatcher ${fusioncatcher_fusions} " : ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fusion_report run $meta.id . $fusionreport_ref $tools --allow-multiple-gene-symbols

    mv fusion_list.tsv ${prefix}.fusionreport.tsv
    mv fusion_list_filtered.tsv ${prefix}.fusionreport_filtered.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusion_report: \$(fusion_report --version)
    END_VERSIONS
    """
}
