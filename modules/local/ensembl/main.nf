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
    tuple val(meta), path("Homo_sapiens.${genome}.${ensembl_version}.gtf")                       , emit: gtf
    tuple val(meta), path("Homo_sapiens.${genome}.${ensembl_version}.dna.primary_assembly.fa")   , emit: fasta
    path "versions.yml"                                                                          , emit: versions


    script:
    """
    if [ ${genome} == 'GRCh37' ]; then
        wget ftp://ftp.ensembl.org/pub/grch37/release-${ensembl_version}/gtf/homo_sapiens/Homo_sapiens.${genome}.87.gtf.gz -O Homo_sapiens.${genome}.${ensembl_version}.gtf.gz
        wget ftp://ftp.ensembl.org/pub/grch37/release-${ensembl_version}/fasta/homo_sapiens/dna/Homo_sapiens.${genome}.dna.primary_assembly.fa.gz -O Homo_sapiens.${genome}.${ensembl_version}.dna.primary_assembly.fa.gz


    else:
        wget ftp://ftp.ensembl.org/pub/release-${ensembl_version}/gtf/homo_sapiens/Homo_sapiens.${genome}.${ensembl_version}.gtf.gz
        wget ftp://ftp.ensembl.org/pub/release-${ensembl_version}/fasta/homo_sapiens/dna/Homo_sapiens.${genome}.dna.primary_assembly.fa.gz -O Homo_sapiens.${genome}.${ensembl_version}.dna.primary_assembly.fa.gz

    gunzip Homo_sapiens.${genome}.${ensembl_version}.gtf.gz
    gunzip Homo_sapiens.${genome}.${ensembl_version}.dna.primary_assembly.fa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3 > versions.yml)
        gunzip: \$(echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')

    END_VERSIONS
    """

    stub:
    """
    touch "Homo_sapiens.${genome}.${ensembl_version}.gtf"
    touch "Homo_sapiens.${genome}.${ensembl_version}.dna.primary_assembly.fa"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3 > versions.yml)
        gunzip: \$(echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

}
