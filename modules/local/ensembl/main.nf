process ENSEMBL_DOWNLOAD {
    tag "ensembl"
    label 'process_low'

    conda "bioconda::gnu-wget=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h5bf99c6_5' :
        'biocontainers/gnu-wget:1.18--h5bf99c6_5' }"

    input:
    val ensembl_version
    val genome
    val meta

    output:
    tuple val(meta), path("Homo_sapiens.${genome}.${ensembl_version}.gtf")                       , emit: gtf
    tuple val(meta), path("Homo_sapiens.${genome}.${ensembl_version}.dna.primary_assembly.fa")   , emit: primary_assembly
    tuple val(meta), path("Homo_sapiens.${genome}.${ensembl_version}.chr.gtf")                   , emit: chrgtf
    path "versions.yml"                                                                          , emit: versions


    script:
    """
    wget ftp://ftp.ensembl.org/pub/release-${ensembl_version}/gtf/homo_sapiens/Homo_sapiens.${params.genome}.${ensembl_version}.gtf.gz
    wget ftp://ftp.ensembl.org/pub/release-${ensembl_version}/fasta/homo_sapiens/dna/Homo_sapiens.${params.genome}.dna.primary_assembly.fa.gz -O Homo_sapiens.${params.genome}.${ensembl_version}.dna.primary_assembly.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-${ensembl_version}/gtf/homo_sapiens/Homo_sapiens.${params.genome}.${ensembl_version}.chr.gtf.gz

    gunzip Homo_sapiens.${params.genome}.${ensembl_version}.gtf.gz
    gunzip Homo_sapiens.${params.genome}.${ensembl_version}.dna.primary_assembly.fa.gz
    gunzip Homo_sapiens.${params.genome}.${ensembl_version}.chr.gtf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3)
        gunzip: \$(echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')

    END_VERSIONS
    """

    stub:
    """
    touch "Homo_sapiens.${genome}.${ensembl_version}.gtf"
    touch "Homo_sapiens.${params.genome}.${ensembl_version}.dna.primary_assembly.fa"
    touch "Homo_sapiens.${params.genome}.${ensembl_version}.chr.gtf"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3)
        gunzip: \$(echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

}
