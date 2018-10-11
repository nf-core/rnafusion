FROM nfcore/rnafusion

ARG FC_GENOME=saccharomyces_cerevisiae

RUN mkdir /data

# Fusioncatcher data
RUN mkdir /data/fc && fusioncatcher-build -g ${FC_GENOME} -o /data/fc

# STAR-Fusion
RUN conda install -c biocore hmmer
RUN mkdir /data/sf
WORKDIR /data/sf
RUN wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz -O Pfam-A.hmm.gz
RUN gunzip Pfam-A.hmm.gz && hmmpress Pfam-A.hmm 
RUN wget https://raw.githubusercontent.com/ndaniel/fusioncatcher/master/bin/get_genome.py -O get_genome.py
RUN sed -i -e 's/import concatenate/from numpy import concatenate/g' get_genome.py
RUN python get_genome.py --organism saccharomyces_cerevisiae --output .
RUN prep_genome_lib.pl --genome_fa genome.fa --gtf /data/fc/organism.gtf --pfarm_db Pfam-A.hmm --CPU 8

# Install nextflow
RUN curl -s https://get.nextflow.io | bash
RUN chmod +x ./nextflow && mv nextflow /usr/local/bin

# Get test data
RUN mkdir /test-data
RUN cd /test-data && wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR612/006/SRR6129616/SRR6129616_1.fastq.gz -O test_1.fastq.gz
RUN cd /test-data && wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR612/006/SRR6129616/SRR6129616_2.fastq.gz -O test_2.fastq.gz

ENV USER travis-ci-test

WORKDIR /
