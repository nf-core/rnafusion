FROM ubuntu:16.04


LABEL authors="rickard.hammaren@scilifelab.se, phil.ewels@scilifelab.se" \
    description="Docker image containing all requirements for NGI-RNAfusion pipeline"

#Fusioncatcher 
RUN apt-get -y update
RUN apt-get -y install \
    build-essential\
    libncurses5-dev \
    default-jdk \
    gawk \
    gcc \
    g++ \
    bzip2 \
    cmake \
    automake \
    gzip \
    zip \
    unzip \
    zlib1g-dev \
    zlib1g \
    wget \
    curl \
    pigz \
    tar \
    parallel \
    libtbb-dev \
    libtbb2 \
    python \
    python-dev \
    python-numpy \
    python-biopython \
    python-xlrd \
    python-openpyxl

RUN wget http://sf.net/projects/fusioncatcher/files/bootstrap.py -O bootstrap.py && python bootstrap.py -t --download -y && rm bootstrap.py



#### STAR-Fusion START 
RUN apt-get update && apt-get install -y gcc g++ perl python automake make \
                                       wget curl libdb-dev \
                       bzip2 zlibc zlib1g zlib1g-dev  default-jre \
                       unzip && \
    apt-get clean

RUN curl -L https://cpanmin.us | perl - App::cpanminus
## set up tool config and deployment area:

ENV SRC /usr/local/src
ENV BIN /usr/local/bin

ENV DATA /usr/local/data
RUN mkdir $DATA


## perl lib installations

RUN cpanm install Set::IntervalTree  # now included w/ STAR-Fusion
RUN cpanm install DB_File
RUN cpanm install URI::Escape
RUN cpanm install Carp::Assert
RUN cpanm install JSON::XS.pm

######################
## Tool installations:
######################



########
# GMAP

RUN GMAP_URL="http://research-pub.gene.com/gmap/src/gmap-gsnap-2017-01-14.tar.gz" && \
    cd $SRC && \
    wget $GMAP_URL && \
    tar xvf gmap-gsnap-2017-01-14.tar.gz && \
    cd gmap-2017-01-14 && ./configure && make && make install


ENV PERL5LIB ${STAR_FUSION_HOME}/PerlLib





###############
## STAR-Fusion:

ENV STAR_FUSION_VERSION=1.1.0

RUN STAR_FUSION_URL="https://github.com/STAR-Fusion/STAR-Fusion/releases/download/v${STAR_FUSION_VERSION}/STAR-Fusion_v${STAR_FUSION_VERSION}.tar.gz" && \
    cd $SRC && \
    wget $STAR_FUSION_URL && \
    tar xvf STAR-Fusion_v${STAR_FUSION_VERSION}.tar.gz && \
    cd STAR-Fusion_v${STAR_FUSION_VERSION} && make

ENV STAR_FUSION_HOME $SRC/STAR-Fusion_v${STAR_FUSION_VERSION}

## FusionInspector now included with STAR-Fusion


########
# Samtools

RUN SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2" && \
   cd $SRC && \
   wget $SAMTOOLS_URL && \
   tar xvf samtools-1.3.1.tar.bz2 && \
   cd samtools-1.3.1/htslib-1.3.1 && ./configure && make && make install && \
   cd ../ && ./configure --without-curses && make && make install


########
# Trinity

ENV TRINITY_VERSION=2.4.0

RUN TRINITY_URL="https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v${TRINITY_VERSION}.tar.gz" && \
   cd $SRC && \
   wget $TRINITY_URL && \
   tar xvf Trinity-v${TRINITY_VERSION}.tar.gz && \
   cd trinityrnaseq-Trinity-v${TRINITY_VERSION} && make


ENV TRINITY_HOME /usr/local/src/trinityrnaseq-Trinity-v${TRINITY_VERSION}

