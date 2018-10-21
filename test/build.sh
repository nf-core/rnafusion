#!/bin/bash

VERSION="nf-core-rnafusion-1.0dev"

# Base environment and tools environment
cp ../environment.yml . && cp -r ../tools .

cat > Dockerfile <<EOF
# This is a testing container with all the tools except fusioncatcher which requires 40GB
# of RAM as a starter. This build usually takes about 1h to finish.
# NOTE: prep_genome_lib.pl ignores CPU parameter so manual edit was applied

FROM nfcore/test-datasets:rnafusion

LABEL authors="rickard.hammaren@scilifelab.se, phil.ewels@scilifelab.se, martin.proks@scilifelab.se" \\
    description="Docker image containing all requirements for nfcore/rnafusion-test container"

ENV USER travis-ci-test

# Get all tools
COPY environment.yml .
RUN conda env create -f /environment.yml python=2.7 && conda clean -a 2> /dev/null
ENV PATH /opt/conda/envs/${VERSION}/bin:\$PATH
EOF

for tool in ./tools/*/; do
    for requirements in $tool*; do
        extension=`echo ${requirements##*.}`
        if [ $extension == "yml" ]
        then
            echo "RUN conda env create -f /tools/${requirements} && conda clean -a 2> /dev/null" >> Dockerfile
        fi
    done
done

cat >> Dockerfile <<EOF

# STAR-Fusion
## perl lib installations
RUN apt-get -y install make build-essential gcc 2> /dev/null
RUN curl -L https://cpanmin.us | perl - App::cpanminus
RUN cpanm install Set::IntervalTree \\
    DB_File::DB_Database \\
    URI::Escape \\
    Set::IntervalTree \\
    Carp::Assert \\
    JSON::XS.pm \\
    PerlIO::gzip
WORKDIR /data/
RUN hmmpress Pfam-A.hmm 
RUN sed -i -e 's/my $CPU = 4;/my $CPU = 2;/' /opt/conda/envs/nf-core-rnafusion-1.0dev/lib/STAR-Fusion/FusionFilter/prep_genome_lib.pl
RUN prep_genome_lib.pl --genome_fa genome.fa --gtf annotation.gtf --pfarm_db Pfam-A.hmm 2> /dev/null

# Install nextflow
WORKDIR /
RUN curl -s https://get.nextflow.io | bash
RUN chmod +x ./nextflow && mv nextflow /usr/local/bin
EOF

# Build image
docker build . -t "rnafusion-test"

# Cleanup
rm ./environment.yml
rm -rf ./tools