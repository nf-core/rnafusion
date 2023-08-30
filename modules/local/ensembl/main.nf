process ENSEMBL_DOWNLOAD {
    tag "ensembl"
    label 'process_low'

    conda "bioconda::gnu-wget=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h5bf99c6_5' :
        'quay.io/biocontainers/gnu-wget:1.18--h5bf99c6_5' }"

    input:
    val ensembl_version
    val genome
    val meta

    output:
    tuple val(meta), path("Homo_sapiens.${genome}.${ensembl_version}.all.fa")        , emit: fasta
    tuple val(meta), path("Homo_sapiens.${genome}.${ensembl_version}.gtf")           , emit: gtf
    tuple val(meta), path("Homo_sapiens.${genome}.${ensembl_version}.chr.gtf")       , emit: chrgtf
    tuple val(meta), path("Homo_sapiens.${genome}.${ensembl_version}.cdna.all.fa.gz"), emit: transcript
    path "versions.yml"                                                                                                                       , emit: versions


    script:
    """
    wget ftp://ftp.ensembl.org/pub/release-${ensembl_version}/fasta/homo_sapiens/dna/Homo_sapiens.${params.genome}.dna.chromosome.{1..22}.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-${ensembl_version}/fasta/homo_sapiens/dna/Homo_sapiens.${params.genome}.dna.chromosome.{MT,X,Y}.fa.gz

    wget ftp://ftp.ensembl.org/pub/release-${ensembl_version}/gtf/homo_sapiens/Homo_sapiens.${params.genome}.${ensembl_version}.gtf.gz
    wget ftp://ftp.ensembl.org/pub/release-${ensembl_version}/gtf/homo_sapiens/Homo_sapiens.${params.genome}.${ensembl_version}.chr.gtf.gz
    wget ftp://ftp.ensembl.org/pub/release-${ensembl_version}/fasta/homo_sapiens/cdna/Homo_sapiens.${params.genome}.cdna.all.fa.gz -O Homo_sapiens.${params.genome}.${ensembl_version}.cdna.all.fa.gz

    gunzip -c Homo_sapiens.${params.genome}.dna.chromosome.* > Homo_sapiens.${params.genome}.${ensembl_version}.all.fa
    gunzip Homo_sapiens.${params.genome}.${ensembl_version}.gtf.gz
    gunzip Homo_sapiens.${params.genome}.${ensembl_version}.chr.gtf.gz



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3 > versions.yml)
    END_VERSIONS
    """

    stub:
    """
    touch "Homo_sapiens.${genome}.${ensembl_version}.all.fa"
    touch "Homo_sapiens.${genome}.${ensembl_version}.gtf"
    touch "Homo_sapiens.${genome}.${ensembl_version}.chr.gtf"
    touch "Homo_sapiens.${genome}.${ensembl_version}.cdna.all.fa.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3 > versions.yml)
    END_VERSIONS
    """

}
