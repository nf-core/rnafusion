process FUSIONCATCHER {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::fusioncatcher=1.33" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "docker.io/clinicalgenomics/fusioncatcher:1.33"
    } else {
        container "docker.io/clinicalgenomics/fusioncatcher:1.33"
    }

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
    """
    fusioncatcher.py \\
        -d $reference \\
        -i $reads \\
        -p $task.cpus \\
        -o . \\
        --skip-blat \\
        $args

    mv final-list_candidate-fusion-genes.txt ${prefix}.fusioncatcher.fusion-genes.txt
    mv summary_candidate_fusions.txt ${prefix}.fusioncatcher.summary.txt
    mv fusioncatcher.log ${prefix}.fusioncatcher.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioncatcher: \$(echo \$(fusioncatcher --version 2>&1)| sed 's/fusioncatcher.py //')
    END_VERSIONS
    touch versions.yml
    """
}
