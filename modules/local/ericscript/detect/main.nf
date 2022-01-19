process ERICSCRIPT {
    tag "eriscript"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publishDir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publishDir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::ericscript=0.5.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ericscript:0.5.5--hdfd78af_5"
    } else {
        container "quay.io/biocontainers/ericscript:0.5.5--hdfd78af_5"
    }

    input:
    tuple val(meta), path(reads)
    path reference

    output:
    path "versions.yml"                                                    , emit: versions
    tuple val(meta), path "./tmp/${prefix}.ericscript.results.filtered.tsv", emit: fusions
    tuple val(meta), path "./tmp/${prefix}.ericscript.results.total.tsv"   , emit: fusions_total

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    ericscript.pl \\
        -db $reference \\
        -name ${prefix}.ericscript \\
        -p $task.cpus \\
        -o ./tmp \\
        $reads \\
        $args

    echo \$(wget -V 2>&1) | grep "GNU Wget" | cut -d" " -f3 > versions.yml
    """
}
