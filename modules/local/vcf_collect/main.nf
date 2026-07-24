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
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('vcf_collect'), val('0.1'), topic: versions, emit: versions_vcf_collect
    tuple val(meta), path("*vcf")    , emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vcf_collect.py --fusioninspector $fusioninspector_tsv --fusionreport $fusionreport_report --fusioninspector_gtf $fusioninspector_gtf_tsv --fusionreport_csv $fusionreport_csv --hgnc $hgnc_ref --sample ${prefix} --out ${prefix}_fusion_data.vcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_fusion_data.vcf
    """
}
