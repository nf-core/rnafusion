FROM openjdk:8

LABEL authors="rickard.hammaren@scilifelab.se, phil.ewels@scilifelab.se" \
    description="Docker image containing all requirements for NGI-RNAfusion pipeline"



#Fusioncatcher 
RUN apt-get -y install \
    build-essential\
    libncurses5-dev \
    default-jdk \
    gawk \
    gcc \
    g++ \
    bzip2 \
    make \
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
RUN wget http://sf.net/projects/fusioncatcher/files/bootstrap.py -O bootstrap.py && python bootstrap.py -t --download -y


#FusionInspector
RUN wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.4.0.tar.gz 
RUN make 
RUN make plugins



