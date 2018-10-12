#!/bin/bash

cp ../environment.yml .
# Get Pfam and test data
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz -O data/Pfam-A.hmm.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR612/006/SRR6129616/SRR6129616_1.fastq.gz -O test-data/test_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR612/006/SRR6129616/SRR6129616_2.fastq.gz -O test-data/test_2.fastq.gz
# Build docker
docker build . -t "rnafusion-test"
