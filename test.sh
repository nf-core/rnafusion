#!/bin/bash

if [ $# -eq 0 ]; then
    echo "No tool provided, testing default pipeline..."
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR612/006/SRR6129616/SRR6129616_1.fastq.gz -O ${TRAVIS_BUILD_DIR}/tests/test_1.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR612/006/SRR6129616/SRR6129616_2.fastq.gz -O ${TRAVIS_BUILD_DIR}/tests/test_2.fastq.gz
    nextflow run ${TRAVIS_BUILD_DIR} -profile test -with-docker nfcore/rnafusion:dev --reads "${TRAVIS_BUILD_DIR}/tests/test_{1,2}.fastq.gz" --genome R64-1-1
  else
    TOOL=$1
    if [ ! -d "./tools/$TOOL" ]; then
        echo "The tool doesn't exist"
        exit 1
    else
        if [ ! -d ${TRAVIS_BUILD_DIR}/tests ]; then
            mkdir ${TRAVIS_BUILD_DIR}/tests
        fi
        cp ${TRAVIS_BUILD_DIR}/tools/$TOOL/* ${TRAVIS_BUILD_DIR}/tests/
        cd ${TRAVIS_BUILD_DIR}/tests && mv Dockerfile.test Dockerfile
        CONTAINER=test-$TOOL
        PROCESS=`echo $TOOL | sed 's/-/_/g'`
        docker build -t $CONTAINER .
        docker run -v ${TRAVIS_BUILD_DIR}:${TRAVIS_BUILD_DIR} $CONTAINER nextflow run ${TRAVIS_BUILD_DIR} -profile test --reads "/test-data/test_{1,2}.fastq.gz" --genome R64-1-1 --$PROCESS --test
    fi
fi