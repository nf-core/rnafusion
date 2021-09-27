// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ERICSCRIPT {
    tag "eriscript"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

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
    path "*.version.txt"                                                   , emit: version
    tuple val(meta), path "./tmp/${prefix}.ericscript.results.filtered.tsv", emit: fusions
    tuple val(meta), path "./tmp/${prefix}.ericscript.results.total.tsv"   , emit: fusions_total
    
    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    ericscript.pl \\
        -db $reference \\
        -name ${prefix}.ericscript \\
        -p $task.cpus \\
        -o ./tmp \\
        $reads \\
        $options.args

    echo \$(wget -V 2>&1) | grep "GNU Wget" | cut -d" " -f3 > ${software}.version.txt
    """
}