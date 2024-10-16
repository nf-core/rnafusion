process GENCODE_DOWNLOAD {
    tag "gencode_download"
    label 'process_low'

    conda "bioconda::gnu-wget=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h5bf99c6_5' :
        'quay.io/biocontainers/gnu-wget:1.18--h5bf99c6_5' }"

    input:
    val genome_gencode_version
    val genome

    output:
    path "*.fa"        , emit: fasta
    path "*.gtf"       , emit: gtf
    path "versions.yml", emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def folder_gencode = genome.contains("38") ? "" : "/${genome}_mapping"
    def gtf_file_name = genome.contains("38") ? "gencode.v${genome_gencode_version}.primary_assembly.annotation.gtf.gz" : "gencode.v${genome_gencode_version}lift${genome_gencode_version}.annotation.gtf.gz"
    """
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${genome_gencode_version}/${folder_gencode}${genome}.primary_assembly.genome.fa.gz -O Homo_sapiens_${genome}_${genome_gencode_version}_dna_primary_assembly.fa.gz
    gunzip Homo_sapiens_${genome}_${genome_gencode_version}_dna_primary_assembly.fa.gz
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${genome_gencode_version}/${folder_gencode}${gtf_file_name} -O Homo_sapiens_${genome}_${genome_gencode_version}.gtf.gz
    gunzip  Homo_sapiens_${genome}_${genome_gencode_version}.gtf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3 > versions.yml)
    END_VERSIONS
    """

    stub:
    """
    touch Homo_sapiens.${genome}.${genome_gencode_version}_dna_primary_assembly.fa
    touch Homo_sapiens.${genome}.${genome_gencode_version}.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo wget -V 2>&1 | grep "GNU Wget" | cut -d" " -f3 > versions.yml)
    END_VERSIONS
    """

}
