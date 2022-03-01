process FUSIONCATCHER {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::fusioncatcher=1.33" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        // container "https://depot.galaxyproject.org/singularity/fusioncatcher:1.30--hdfd78af_1"
        container "quay.io/biocontainers/fusioncatcher:1.33--hdfd78af_2"
    } else {
        container "https://depot.galaxyproject.org/singularity/fusioncatcher:1.30--0"
        // container "quay.io/biocontainers/fusioncatcher:1.33--hdfd78af_1"
    }

    input:
    tuple val(meta), path(fasta)
    path reference

    output:
    tuple val(meta), path("*.fusioncatcher.fusion-genes.txt")   , optional:true, emit: fusions
    tuple val(meta), path("*.fusioncatcher.summary.txt")        , optional:true, emit: summary
    tuple val(meta), path("*.fusioncatcher.log")                               , emit: log
    path "versions.yml"                                                        , emit: versions

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
    touch versions.yml
    """
}


    // mv final-list_candidate-fusion-genes.txt ${prefix}.fusioncatcher.fusion-genes.txt
    // mv summary_candidate_fusions.txt ${prefix}.fusioncatcher.summary.txt
    // mv fusioncatcher.log ${prefix}.fusioncatcher.log
    //  fusioncatcher --version | sed 's/fusioncatcher.py //' > versions.yml
    // touch ${prefix}.fusioncatcher.fusion-genes.txt
    // touch ${prefix}.fusioncatcher.summary.txt
    // touch ${prefix}.fusioncatcher.log
