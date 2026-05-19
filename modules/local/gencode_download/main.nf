process GENCODE_DOWNLOAD {
    tag "gencode_download"
    label 'process_low'

    conda "bioconda::gnu-wget=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3b/3b54fa9135194c72a18d00db6b399c03248103f87e43ca75e4b50d61179994b3/data' :
        'community.wave.seqera.io/library/wget:1.21.4--8b0fcde81c17be5e' }"

    input:
    val genome_gencode_version
    val genome

    output:
    path "*.fa"        , emit: fasta
    path "*.gtf"       , emit: gtf
    tuple val("${task.process}"), val('wget'), eval("wget --version | head -1 | cut -d ' ' -f 3"), topic: versions, emit: versions_wget


    when:
    task.ext.when == null || task.ext.when

    script:
    def folder_gencode = genome.contains("38") ? "" : "/${genome}_mapping"
    def gtf_file_name = genome.contains("38") ? "gencode.v${genome_gencode_version}.primary_assembly.annotation.gtf.gz" : "gencode.v${genome_gencode_version}lift${genome_gencode_version}.annotation.gtf.gz"
    """
    wget --no-check-certificate ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${genome_gencode_version}/${folder_gencode}${genome}.primary_assembly.genome.fa.gz -O Homo_sapiens_${genome}_${genome_gencode_version}_dna_primary_assembly.fa.gz
    gunzip Homo_sapiens_${genome}_${genome_gencode_version}_dna_primary_assembly.fa.gz
    wget --no-check-certificate ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${genome_gencode_version}/${folder_gencode}${gtf_file_name} -O Homo_sapiens_${genome}_${genome_gencode_version}.gtf.gz
    gunzip  Homo_sapiens_${genome}_${genome_gencode_version}.gtf.gz
    """

    stub:
    """
    touch Homo_sapiens.${genome}.${genome_gencode_version}_dna_primary_assembly.fa
    touch Homo_sapiens.${genome}.${genome_gencode_version}.gtf
    """

}
