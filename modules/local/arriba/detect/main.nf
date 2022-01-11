process ARRIBA {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publishDir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publishDir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::arriba=2.1.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/arriba:2.1.0--h3198e80_1"
    } else {
        container "quay.io/biocontainers/arriba:2.1.0--h3198e80_1"
    }

    input:
    tuple val(meta), path(reads)
    path fasta
    path gtf
    path index

    output:
    tuple val(meta), path("*.fusions.tsv")          , emit: fusions
    tuple val(meta), path("*.fusions.discarded.tsv"), emit: fusions_fail
    path "versions.yml"                             , emit: versions

    script:
    def prefix    = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    STAR \\
        --genomeDir $index \\
        --readFilesIn $reads \\
        --runThreadN $task.cpus \\
        $args |

    tee Aligned.out.bam |

    arriba \\
        -x /dev/stdin \\
        -a $fasta \\
        -g $gtf \\
        -o ${prefix}.arriba.fusions.tsv \\
        -O ${prefix}.arriba.fusions.discarded.tsv \\
        $args2

    echo \$(arriba -h | grep 'Version:' 2>&1) |  sed 's/Version:\s//' > versions.yml
    """
}
