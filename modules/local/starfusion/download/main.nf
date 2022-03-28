process STARFUSION_DOWNLOAD {
    tag 'star-fusion'

    conda (params.enable_conda ? "bioconda::dfam=3.3 bioconda::hmmer=3.3.2 bioconda::star-fusion=1.10.0 bioconda::trinity=date.2011_11_2 bioconda::samtools=1.9 bioconda::star=2.7.8a" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "docker.io/trinityctat/starfusion:1.10.1"
    } else {
        container "docker.io/trinityctat/starfusion:1.10.1"
    }

    input:
    path fasta
    path gtf

    output:
    path "*"  , emit: reference

    script:
    """
    wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz --no-check-certificate

    tar xvf CTAT_resource_lib.tar.gz

    rm CTAT_resource_lib.tar.gz
    """
}
