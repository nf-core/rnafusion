process GFFREAD {
    tag "$gff"
    label 'process_low'

    conda "bioconda::gffread=0.12.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.1--h8b12597_0' :
        'biocontainers/gffread:0.12.1--h8b12597_0' }"

    input:
    tuple val(meta), path(gff)
    tuple val(meta2), path(fasta), path(fai)

    output:
    path "*_gffread_output.gtf", emit: gtf,        optional: true
    path "*_gffread_output.fa"  , emit: tr_fasta,   optional: true
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def genome = fasta ? "-g ${fasta}" : ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outp   = args.contains("-w") ? "-w ${prefix}_gffread_output.fa" : "-o ${prefix}_gffread_output.gtf"

    """
    gffread \\
        $gff \\
        $genome \\
        $outp
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_gffread"
    """
    touch ${prefix}_gffread_output.gtf
    touch ${prefix}_gffread_output.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """
}
