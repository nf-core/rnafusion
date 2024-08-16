process FUSIONCATCHER {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::fusioncatcher=1.33"
    container "nf-core/fusioncatcher:1.33"

    input:
    tuple val(meta), path(fasta)
    path reference

    output:
    tuple val(meta), path("*.fusioncatcher.fusion-genes.txt")   , optional:true  , emit: fusions
    tuple val(meta), path("*.fusioncatcher.summary.txt")        , optional:true  , emit: summary
    tuple val(meta), path("*.fusioncatcher.log")                                 , emit: log
    path "versions.yml"                                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads = fasta.toString().replace(" ", ",")
    def single_end = meta.single_end ? "--single-end" : ""
    """
    fusioncatcher.py \\
        -d $reference \\
        -i $reads \\
        -p $task.cpus \\
        -o . \\
        --skip-blat \\
        $single_end \\
        $args

    mv final-list_candidate-fusion-genes.txt ${prefix}.fusioncatcher.fusion-genes.txt
    mv summary_candidate_fusions.txt ${prefix}.fusioncatcher.summary.txt
    mv fusioncatcher.log ${prefix}.fusioncatcher.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioncatcher: \$(echo \$(fusioncatcher.py --version 2>&1)| sed 's/fusioncatcher.py //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.fusioncatcher.fusion-genes.txt
    touch ${prefix}.fusioncatcher.summary.txt
    touch ${prefix}.fusioncatcher.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioncatcher: \$(echo \$(fusioncatcher.py --version 2>&1)| sed 's/fusioncatcher.py //')
    END_VERSIONS
    """
}
