process VCF_COLLECT {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::gtfparse=2.0.1"
    container "quay.io/biocontainers/gtfparse:2.0.1--pyh7cba7a3_1"

    input:
    tuple val(meta), path(fusioninspector_tsv), path(fusioninspector_gtf_tsv), path(fusionreport_report)
    tuple val(meta2),  path(hgnc_ref)
    tuple val(meta3),  path(hgnc_date)

    output:
    path "versions.yml"              , emit: versions
    tuple val(meta), path("*vcf")    , emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vcf_collect.py --fusioninspector $fusioninspector_tsv --fusionreport $fusionreport_report --fusioninspector_gtf $fusioninspector_gtf_tsv --hgnc $hgnc_ref --sample ${prefix} --out ${prefix}_fusion_data.vcf

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
