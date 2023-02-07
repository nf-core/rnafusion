process STARFUSION_BUILD {
    tag 'star-fusion'

    conda "bioconda::dfam=3.3 bioconda::hmmer=3.3.2 bioconda::star-fusion=1.10.0 bioconda::trinity=date.2011_11_2 bioconda::samtools=1.9 bioconda::star=2.7.8a"
    container "docker.io/trinityctat/starfusion:1.10.1"

    input:
    path fasta
    path gtf

    output:
    path "*"  , emit: reference

    script:
    def binPath = ( workflow.conda.enabled = true ? "prep_genome_lib.pl" : "/usr/local/src/STAR-Fusion/ctat-genome-lib-builder/prep_genome_lib.pl" )
    """
    export TMPDIR=/tmp
    wget http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam34.0/Pfam-A.hmm.gz --no-check-certificate
    wget https://github.com/FusionAnnotator/CTAT_HumanFusionLib/releases/download/v0.3.0/fusion_lib.Mar2021.dat.gz -O CTAT_HumanFusionLib_Mar2021.dat.gz --no-check-certificate
    wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/AnnotFilterRule.pm -O AnnotFilterRule.pm --no-check-certificate
    wget https://www.dfam.org/releases/Dfam_3.4/infrastructure/dfamscan/homo_sapiens_dfam.hmm --no-check-certificate
    wget https://www.dfam.org/releases/Dfam_3.4/infrastructure/dfamscan/homo_sapiens_dfam.hmm.h3f --no-check-certificate
    wget https://www.dfam.org/releases/Dfam_3.4/infrastructure/dfamscan/homo_sapiens_dfam.hmm.h3i --no-check-certificate
    wget https://www.dfam.org/releases/Dfam_3.4/infrastructure/dfamscan/homo_sapiens_dfam.hmm.h3m --no-check-certificate
    wget https://www.dfam.org/releases/Dfam_3.4/infrastructure/dfamscan/homo_sapiens_dfam.hmm.h3p --no-check-certificate
    gunzip Pfam-A.hmm.gz && hmmpress Pfam-A.hmm
    $binPath \\
        --genome_fa $fasta \\
        --gtf $gtf \\
        --annot_filter_rule AnnotFilterRule.pm \\
        --fusion_annot_lib CTAT_HumanFusionLib_Mar2021.dat.gz \\
        --pfam_db Pfam-A.hmm \\
        --dfam_db homo_sapiens_dfam.hmm \\
        --max_readlength $params.read_length \\
        --CPU $task.cpus
    """

    stub:
    touch

}
