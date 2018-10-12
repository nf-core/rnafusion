#!/bin/bash

cp ../environment.yml .
mkdir data/
wget ftp://ftp.ensembl.org/pub/current_fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz -O data/genome.fa.gz
wget ftp://ftp.ensembl.org/pub/release-94/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.94.gtf.gz -O data/annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz -O data/Pfam-A.hmm.gz
docker build . -t "rnafusion-test"
