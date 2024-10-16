process FUSIONCATCHER_BUILD {
    tag "fusioncatcher_build"
    label 'process_medium'

    conda "bioconda::fusioncatcher=1.33"
    container "docker.io/clinicalgenomics/fusioncatcher:1.33"

    input:
    val genome_gencode_version

    output:
    path "human_v${genome_gencode_version}"  , emit: reference
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''
    """
    fusioncatcher-build.py \\
        -g homo_sapiens \\
        -o human_v${genome_gencode_version} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioncatcher: \$(echo \$(fusioncatcher --version 2>&1))
    END_VERSIONS
    """

    stub:
    """
    mkdir human_v${genome_gencode_version}
    touch human_v${genome_gencode_version}/ensembl_fully_overlapping_genes.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioncatcher: \$(echo \$(fusioncatcher --version 2>&1))
    END_VERSIONS
    """
}
