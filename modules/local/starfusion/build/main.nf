process STARFUSION_BUILD {
    tag 'star-fusion'

    conda "bioconda::dfam=3.7 bioconda::hmmer=3.4 bioconda::star-fusion=1.13.0 bioconda::trinity=2.15.1 bioconda::samtools=1.19.2 bioconda::star=2.7.11b"
    container "docker.io/trinityctat/starfusion:1.13.0"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(gtf)

    output:
    path "*"  , emit: reference

    script:
    def binPath = (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1)  ? "prep_genome_lib.pl" : "/usr/local/src/STAR-Fusion/ctat-genome-lib-builder/prep_genome_lib.pl"
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR-Fusion: \$(STAR-Fusion --version 2>&1 | grep -i 'version' | sed 's/STAR-Fusion version: //')
    END_VERSIONS
    """

    stub:
    """
    mkdir ctat_genome_lib_build_dir
    touch ref_annot.cdna.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR-Fusion: \$(STAR-Fusion --version 2>&1 | grep -i 'version' | sed 's/STAR-Fusion version: //')
    END_VERSIONS
    """

}
