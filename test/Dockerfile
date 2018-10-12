# Species: saccharomyces_cerevisiae
# Dataset taken from: https://www.ebi.ac.uk/ena/data/view/SRX3242046

FROM nfcore/base

LABEL authors="rickard.hammaren@scilifelab.se, phil.ewels@scilifelab.se, martin.proks@scilifelab.se" \
    description="Docker image containing all requirements for NGI-RNAfusion pipeline"

ENV USER travis-ci-test

# Get all tools
COPY environment.yml .
RUN conda env create -f /environment.yml python=2.7 && conda clean -a
ENV PATH /opt/conda/envs/nf-core-rnafusion-1.0dev/bin:$PATH

# Data directory for genome and annotation
COPY ./data /data
# Get test data
COPY ./test-data /test-data

WORKDIR /data
# Genome and GTF
RUN gzip -d -f -c genome.fa.gz > 'genome.fa'
RUN gunzip annotation.gtf.gz
# STAR-Fusion
RUN gunzip Pfam-A.hmm.gz && hmmpress Pfam-A.hmm 
RUN prep_genome_lib.pl --genome_fa genome.fa --gtf annotation.gtf --pfarm_db Pfam-A.hmm --CPU 2 >/dev/null 2>&1

# Install nextflow
WORKDIR /
RUN curl -s https://get.nextflow.io | bash
RUN chmod +x ./nextflow && mv nextflow /usr/local/bin
