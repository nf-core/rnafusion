# Species: saccharomyces_cerevisiae
# Dataset taken from: https://www.ebi.ac.uk/ena/data/view/SRX3242046

FROM nfcore/rnafusion

ARG FC_GENOME=saccharomyces_cerevisiae
ENV USER travis-ci-test

# Get Genome
RUN mkdir /data
WORKDIR /data
RUN wget https://raw.githubusercontent.com/ndaniel/fusioncatcher/master/bin/get_genome.py -O get_genome.py
RUN sed -i -e 's/import concatenate/from numpy import concatenate/g' get_genome.py
RUN python get_genome.py --organism ${FC_GENOME} --output .

# Get GTF
RUN wget ftp://ftp.ensembl.org/pub/release-94/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.94.gtf.gz -O annotation.gtf.gz
RUN gunzip annotation.gtf.gz

# STAR-Fusion
RUN mkdir /data/sf
WORKDIR /data/sf

RUN conda install -c biocore hmmer
RUN wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz -O Pfam-A.hmm.gz
RUN gunzip Pfam-A.hmm.gz && hmmpress Pfam-A.hmm 
RUN prep_genome_lib.pl --genome_fa /data/genome.fa --gtf /data/annotation.gtf --pfarm_db Pfam-A.hmm --CPU 8

# Install nextflow
WORKDIR /
RUN conda install -c anaconda graphviz 
RUN curl -s https://get.nextflow.io | bash
RUN chmod +x ./nextflow && mv nextflow /usr/local/bin

# Get test data
RUN mkdir /test-data
RUN cd /test-data && wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR612/006/SRR6129616/SRR6129616_1.fastq.gz -O test_1.fastq.gz
RUN cd /test-data && wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR612/006/SRR6129616/SRR6129616_2.fastq.gz -O test_2.fastq.gz
