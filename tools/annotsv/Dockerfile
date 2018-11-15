FROM nfcore/base

LABEL authors="martin.proks@scilifelab.se" description="Docker image for AnnotSV tool"

ENV VERSION="1.2"
RUN apt-get -y update
RUN apt-get -y install bedtools tar tcl

WORKDIR /usr/local/src
RUN wget http://lbgi.fr/AnnotSV/Sources/AnnotSV_${VERSION}.tar.gz -O AnnotSV_${VERSION}.tar.gz && \
    tar -xvzf AnnotSV_${VERSION}.tar.gz

ENV ANNOTSV /usr/local/src/AnnotSV_${VERSION}/
ENV PATH /usr/local/src/AnnotSV_${VERSION}/bin:$PATH
WORKDIR /
