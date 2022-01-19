process STARFUSION_DOWNLOAD {
    tag 'star-fusion'
    label 'process_high'
    publishDir "${params.genomes_base}",
        mode: params.publishDir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }

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
    def binPath = if (workflow.containerEngine == "singularity" || workflow.containerEngine == "docker") ? "/usr/local/src/STAR-Fusion/ctat-genome-lib-builder/" : ""
    """
    export TMPDIR=/tmp

    wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam34.0/Pfam-A.hmm.gz
    wget https://github.com/FusionAnnotator/CTAT_HumanFusionLib/releases/download/v0.3.0/fusion_lib.Mar2021.dat.gz -O CTAT_HumanFusionLib_Mar2021.dat.gz
    wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/AnnotFilterRule.pm -O AnnotFilterRule.pm
    wget https://www.dfam.org/releases/Dfam_3.4/infrastructure/dfamscan/homo_sapiens_dfam.hmm
    wget https://www.dfam.org/releases/Dfam_3.4/infrastructure/dfamscan/homo_sapiens_dfam.hmm.h3f
    wget https://www.dfam.org/releases/Dfam_3.4/infrastructure/dfamscan/homo_sapiens_dfam.hmm.h3i
    wget https://www.dfam.org/releases/Dfam_3.4/infrastructure/dfamscan/homo_sapiens_dfam.hmm.h3m
    wget https://www.dfam.org/releases/Dfam_3.4/infrastructure/dfamscan/homo_sapiens_dfam.hmm.h3p

    gunzip Pfam-A.hmm.gz && hmmpress Pfam-A.hmm

    ${binPath}prep_genome.pl \\
        --genome_fa $fasta \\
        --gtf $gtf \\
        --annot_filter_rule AnnotFilterRule.pm \\
        --fusion_annot_lib CTAT_HumanFusionLib_Mar2021.dat.gz \\
        --pfam_db Pfam-A.hmm \\
        --dfam_db homo_sapiens_dfam.hmm \\
        --max_readlength $params.read_length \\
        --CPU $task.cpus
    """

}
