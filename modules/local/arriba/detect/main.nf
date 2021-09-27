// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ARRIBA {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

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
    path "*.version.txt"                            , emit: version

    script:
    def software  = getSoftwareName(task.process)
    def prefix    = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    STAR \\
	    --genomeDir $index \\
	    --readFilesIn $reads \\
	    --runThreadN $task.cpus \\
        $options.args |

    tee Aligned.out.bam |

    arriba \\
        -x /dev/stdin \\
        -a $fasta \\
        -g $gtf \\
        -o ${prefix}.arriba.fusions.tsv \\
        -O ${prefix}.arriba.fusions.discarded.tsv \\
        $options.args2

    echo \$(arriba -h | grep 'Version:' 2>&1) |  sed 's/Version:\s//' > ${software}.version.txt
    """
}
