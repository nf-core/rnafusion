FROM nfcore/base

LABEL authors="rickard.hammaren@scilifelab.se, phil.ewels@scilifelab.se, martin.proks@scilifelab.se" \
    description="Docker image containing all requirements for NGI-RNAfusion pipeline"

LABEL VERSION="1.00 (045a8af)"
RUN apt-get -y update
RUN apt-get -y install \
    automake \
    build-essential \
    bzip2 \
    cmake \
    curl \
    default-jdk \
    default-jre \
    g++ \
    gawk \
    gcc \
    gzip \
    libc6-dev \
    libdb-dev \
    libncurses5-dev \
    libtbb-dev \
    libtbb2 \
    make \
    parallel \
    pigz \
    python \
    python-biopython \
    python-dev \
    python-numpy \
    python-openpyxl \
    python-xlrd \
    tar \
    unzip \
    wget \
    zip \
    zlib1g \
    zlib1g-dev \
    zlibc

## perl lib installations
RUN curl -L https://cpanmin.us | perl - App::cpanminus
RUN cpanm install Set::IntervalTree \
    DB_File \
    URI::Escape\
    Set::IntervalTree \
    Carp::Assert \
    JSON::XS.pm \
    PerlIO::gzip

WORKDIR /usr/local/src

RUN pip install numpy biopython xlrd openpyxl
RUN mkdir fusioncatcher && cd fusioncatcher \
    && wget https://raw.githubusercontent.com/ndaniel/fusioncatcher/045a8af2a8b8e8f3fddd6bfcaf1bd1064cf2cc15/bin/bootstrap.py -O bootstrap.py \
    && wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat -O blat \
    && python bootstrap.py -t -y \
    && chmod +x blat \
    && rm bootstrap.py

WORKDIR /
ENV PATH /usr/local/src/fusioncatcher/:$PATH
