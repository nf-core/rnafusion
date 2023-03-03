process FUSIONREPORT {
    tag "$meta.id"
    label 'process_medium'

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda "bioconda::star=2.7.9a"
    container "docker.io/clinicalgenomics/fusion-report:2.1.5p2"


    input:
    tuple val(meta), path(reads), path(arriba_fusions), path(pizzly_fusions), path(squid_fusions), path(starfusion_fusions),  path(fusioncatcher_fusions)
    path(fusionreport_ref)

    output:
    path "versions.yml"                                                 , emit: versions
    tuple val(meta), path("*fusionreport.tsv")                          , emit: fusion_list
    tuple val(meta), path("*fusionreport_filtered.tsv")                 , emit: fusion_list_filtered
    tuple val(meta), path("*.html")                                     , emit: report
    tuple val(meta), path("*.csv")                       , optional:true, emit: csv
    tuple val(meta), path("*.json")                      , optional:true, emit: json

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def tools = params.arriba || params.all         ? "--arriba ${arriba_fusions} " : ''
    tools    += params.pizzly || params.all         ? "--pizzly ${pizzly_fusions} " : ''
    tools    += params.squid  || params.all         ? "--squid ${squid_fusions} " : ''
    tools    += params.starfusion  || params.all    ? "--starfusion ${starfusion_fusions} " : ''
    tools    += params.fusioncatcher  || params.all ? "--fusioncatcher ${fusioncatcher_fusions} " : ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fusion_report run $meta.id . $fusionreport_ref $tools --allow-multiple-gene-symbols $args $args2

    mv fusion_list.tsv ${prefix}.fusionreport.tsv
    mv fusion_list_filtered.tsv ${prefix}.fusionreport_filtered.tsv
    mv fusions.csv ${prefix}.fusions.csv
    mv fusions.json ${prefix}.fusions.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusion_report: \$(fusion_report --version | sed 's/fusion-report //')
        fusion_report DB retrieval: \$(cat $fusionreport_ref/DB-timestamp.txt)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fusionreport_filtered.tsv
    touch ${prefix}.fusionreport.tsv
    touch index.html
    touch ${prefix}.fusions.csv
    touch ${prefix}.fusions.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusion_report: \$(fusion_report --version | sed 's/fusion-report //')
    END_VERSIONS
    """
}
