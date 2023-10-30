process VCF_COLLECT {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda:: gtfparse =2.0.1"
    container "quay.io/biocontainers/gtfparse:2.0.1--pyh7cba7a3_0"

    input:
    tuple val(meta), path(tsv), path(out_gtf), path(report)
    path hgnc_ref
    path hgnc_date

    output:
    path "versions.yml"              , emit: versions
    tuple val(meta), path("*vcf")    , emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vcf_collect.py --fusioninspector $tsv --fusionreport $report --fusioninspector_gtf $out_gtf --hgnc $hgnc_ref --sample ${prefix} --out ${prefix}.vcf

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
