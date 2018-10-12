# Species: saccharomyces_cerevisiae
# Dataset taken from: https://www.ebi.ac.uk/ena/data/view/SRX3242046

FROM nfcore/base

LABEL authors="rickard.hammaren@scilifelab.se, phil.ewels@scilifelab.se, martin.proks@scilifelab.se" \
    description="Docker image containing all requirements for NGI-RNAfusion pipeline"

ARG FC_GENOME=saccharomyces_cerevisiae
ENV USER travis-ci-test

# Get all tools
COPY environment.yml .
RUN conda env create -f /environment.yml python=2.7 && conda clean -a
RUN conda install -c biocore hmmer pyyaml numpy graphviz
ENV PATH /opt/conda/envs/nf-core-rnafusion-1.0dev/bin:$PATH

# Get Genome
RUN mkdir /data
WORKDIR /data
RUN wget ftp://ftp.ensembl.org/pub/current_fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
RUN gzip -d -f -c './Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz' > 'genome.fa'

# Get GTF
RUN wget ftp://ftp.ensembl.org/pub/release-94/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.94.gtf.gz -O annotation.gtf.gz
RUN gunzip annotation.gtf.gz

# STAR-Fusion
RUN mkdir /data/sf
WORKDIR /data/sf

RUN wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz -O Pfam-A.hmm.gz
RUN gunzip Pfam-A.hmm.gz && hmmpress Pfam-A.hmm 
RUN prep_genome_lib.pl --genome_fa /data/genome.fa --gtf /data/annotation.gtf --pfarm_db Pfam-A.hmm --CPU 8

# Install nextflow
WORKDIR /
RUN curl -s https://get.nextflow.io | bash
RUN chmod +x ./nextflow && mv nextflow /usr/local/bin

# Get test data
RUN mkdir /test-data
RUN cd /test-data && wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR612/006/SRR6129616/SRR6129616_1.fastq.gz -O test_1.fastq.gz
RUN cd /test-data && wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR612/006/SRR6129616/SRR6129616_2.fastq.gz -O test_2.fastq.gz
