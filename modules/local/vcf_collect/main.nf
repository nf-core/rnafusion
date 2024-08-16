process VCF_COLLECT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(fusioninspector_tsv), path(fusioninspector_gtf_tsv), path(fusionreport_report), path(fusionreport_csv)
    tuple val(meta2),  path(hgnc_ref)
    tuple val(meta3),  path(hgnc_date)

    output:
    path "versions.yml"              , emit: versions
    tuple val(meta), path("*vcf.gz") , emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vcf_collect.py --fusioninspector $fusioninspector_tsv --fusionreport $fusionreport_report --fusioninspector_gtf $fusioninspector_gtf_tsv --fusionreport_csv $fusionreport_csv --hgnc $hgnc_ref --sample ${prefix} --out ${prefix}_fusion_data.vcf
    gzip ${prefix}_fusion_data.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        HGNC DB retrieval: \$(cat $hgnc_date)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
