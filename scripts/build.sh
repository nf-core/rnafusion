#!/bin/bash

PREFIX="nfcore/rnafusion"

create_container() {
    TOOL_PATH=$1
    VERSION="$(cat $TOOL_PATH/environment.yml | grep "name:" | cut -d":" -f2 | cut -d "_" -f2)"
    TOOL_NAME=`basename $TOOL_PATH`
    CONTAINER_NAME="${PREFIX}:${TOOL_NAME}_${VERSION}"
    echo "Building [$CONTAINER_NAME]"
    docker build $TOOL_PATH -t $CONTAINER_NAME
    docker push $CONTAINER_NAME
}

if [ $# -eq 0 ]; then
    echo "No tool name specified!"
    echo "Run scripts/build.sh -h for help"
    exit 1
fi

if [ $1 == "-h" ]; then
    echo "Utility for building docker containers from tools/"
    echo "Usage: scripts/build.sh [options]"
    echo 
    echo "Options:"
    echo "  all             build all tools including main image"
    echo "  <tool name>     builds specific tool"
    echo "                  Example: sh scripts/build.sh ericscript"
    exit 0
fi

if [ $1 == "all" ]; then
    for TOOL in containers/*/; do
        create_container `pwd`/$TOOL ${TOOL%?}
    done
    # Build main container
    VERSION="$(cat nextflow.config | grep -m1 "container" | cut -d":" -f2 | cut -d "'" -f1)"
    CONTAINER_NAME=$PREFIX:$VERSION
    echo "Building [$CONTAINER_NAME]"
    docker build . -t $CONTAINER_NAME
    docker push $CONTAINER_NAME
else
    TOOL=$1
    TOOL_PATH="$(pwd)/containers/$TOOL"
    if [ ! -d $TOOL_PATH ]; then
        echo "The tool doesn't exist"
        exit 1
    else
        create_container $TOOL_PATH
    fi
fi
