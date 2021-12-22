process STAR_GENOMEGENERATE {
    tag "$fasta"
    label 'process_high'
    publishDir "${params.genomes_base}",
        mode: params.publishDir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }

    conda (params.enable_conda ? "bioconda::star=2.7.8a" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/star:2.7.8a--h9ee0642_1"
    } else {
        container "quay.io/biocontainers/star:2.7.8a--h9ee0642_1"
    }

    input:
    path fasta
    path gtf

    output:
    path "star/*"         , emit: index
    path "versions.yml"   , emit: versionss

    script:
    def memory   = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    def args     = options.args.tokenize()
    """
    mkdir star
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star/ \\
        --genomeFastaFiles $fasta \\
        --sjdbGTFfile $gtf \\
        --sjdbOverhang ${params.read_length - 1} \\
        --runThreadN $task.cpus \\
        $memory \\
        $options.args

    STAR --version | sed -e "s/STAR_//g" > versions.yml
    """
}
