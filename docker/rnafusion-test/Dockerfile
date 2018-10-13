# This is a testing container with all the tools except fusioncatcher which requires 40GB
# of RAM as a starter. This build usually takes about 1h to finish.
# NOTE: prep_genome_lib.pl ignores CPU parameter so manual edit was applied

FROM nfcore/rnafusion-test-dataset

LABEL authors="rickard.hammaren@scilifelab.se, phil.ewels@scilifelab.se, martin.proks@scilifelab.se" \
    description="Docker image containing all requirements for nfcore/rnafusion-test container"

ENV USER travis-ci-test

# Get all tools
COPY environment.yml .
RUN conda env create -f /environment.yml python=2.7 && conda clean -a 2> /dev/null
ENV PATH /opt/conda/envs/nf-core-rnafusion-1.0dev/bin:$PATH

# STAR-Fusion
## perl lib installations
RUN apt-get -y install make build-essential gcc 2> /dev/null
RUN curl -L https://cpanmin.us | perl - App::cpanminus
RUN cpanm install Set::IntervalTree \
    DB_File::DB_Database \
    URI::Escape\
    Set::IntervalTree \
    Carp::Assert \
    JSON::XS.pm \
    PerlIO::gzip
WORKDIR /data/
RUN hmmpress Pfam-A.hmm 
RUN sed -i -e 's/my $CPU = 4;/my $CPU = 2;/' /opt/conda/envs/nf-core-rnafusion-1.0dev/lib/STAR-Fusion/FusionFilter/prep_genome_lib.pl
RUN prep_genome_lib.pl --genome_fa genome.fa --gtf annotation.gtf --pfarm_db Pfam-A.hmm 2> /dev/null

# Install nextflow
WORKDIR /
RUN curl -s https://get.nextflow.io | bash
RUN chmod +x ./nextflow && mv nextflow /usr/local/bin
