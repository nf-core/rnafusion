FROM nfcore/base

LABEL authors="rickard.hammaren@scilifelab.se, phil.ewels@scilifelab.se, martin.proks@scilifelab.se" \
    description="Fake volume for nfcore/rnafusion"

# Download data
RUN mkdir /data && mkdir /test-data

WORKDIR /data
RUN wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz -O Pfam-A.hmm.gz
RUN wget ftp://ftp.ensembl.org/pub/current_fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz -O genome.fa.gz
RUN wget ftp://ftp.ensembl.org/pub/release-94/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.94.gtf.gz -O annotation.gtf.gz

RUN gzip -d -f -c genome.fa.gz > genome.fa
RUN gunzip annotation.gtf.gz
RUN gunzip Pfam-A.hmm.gz

# Download test-data
RUN wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR612/006/SRR6129616/SRR6129616_1.fastq.gz -O /test-data/test_1.fastq.gz
RUN wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR612/006/SRR6129616/SRR6129616_2.fastq.gz -O /test-data/test_2.fastq.gz

WORKDIR /
