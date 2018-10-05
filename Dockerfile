FROM ubuntu:16.04

LABEL authors="rickard.hammaren@scilifelab.se, phil.ewels@scilifelab.se, martin.proks@scilifelab.se" \
    description="Docker image containing all requirements for NGI-RNAfusion pipeline"

#Fusioncatcher 
RUN apt-get -y update
RUN apt-get -y install \
    build-essential\
    libncurses5-dev \
    default-jdk \
    default-jre \
    gawk \
    gcc \
    g++ \
    bzip2 \
    cmake \
    automake \
    make \
    gzip \
    zip \
    unzip \
    zlibc \
    zlib1g-dev \
    zlib1g \
    wget \
    curl \
    pigz \
    tar \
    parallel \
    libtbb-dev \
    libtbb2 \
    libdb-dev \
    python \
    python-dev \
    python-numpy \
    python-biopython \
    python-xlrd \
    python-openpyxl 

ENV SRC /usr/local/src
WORKDIR $SRC

## perl lib installations
RUN curl -L https://cpanmin.us | perl - App::cpanminus
RUN cpanm install Set::IntervalTree \
    DB_File \
    URI::Escape\
    Set::IntervalTree \
    Carp::Assert \
    JSON::XS.pm \
    PerlIO::gzip

# Conda
RUN wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda2.sh && \
    sh miniconda2.sh -b -p $SRC/miniconda2
ENV PATH=$SRC/miniconda2/bin/:$PATH

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN /bin/bash -c "source ${SRC}/miniconda2/bin/activate"
ENV PATH /usr/local/src/miniconda2/envs/nf-core-rnafusion-1.0dev/bin:$PATH

# Fusion Catcher
RUN mkdir fusioncatcher && cd fusioncatcher && \
    wget http://sf.net/projects/fusioncatcher/files/bootstrap.py -O bootstrap.py && \
    python bootstrap.py -t -y && \
    rm bootstrap.py
ENV PATH=$SRC/fusioncatcher/:${PATH}

WORKDIR /