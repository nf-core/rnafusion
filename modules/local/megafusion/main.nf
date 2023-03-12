process MEGAFUSION {
    tag "$meta.id"
    label 'process_single'

    // Note: 2.7X indices incompatible with AWS iGenomes.
    container "docker.io/clinicalgenomics/megafusion:1.0.0"


    input:
    tuple val(meta), path(tsv)

    output:
    path "versions.yml"              , emit: versions
    tuple val(meta), path("*vcf")    , emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python MegaFusion.py --json json/FusionInspector.json --fusion $tsv > ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MegaFusion: 1.0.0
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MegaFusion: 1.0.0
    END_VERSIONS
    """
}
