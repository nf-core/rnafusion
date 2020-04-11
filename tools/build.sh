#!/bin/bash

PREFIX="nfcore/rnafusion"

create_container() {
    TOOL_PATH=$1
    VERSION="$(cat $TOOL_PATH/environment.yml | grep "name:" | cut -d":" -f2 | cut -d " " -f2)"
    CONTAINER_NAME="$PREFIX:$VERSION"
    echo "Building [$CONTAINER_NAME]"
    docker build $TOOL_PATH -t $CONTAINER_NAME
    docker push $CONTAINER_NAME
}

if [ $# -eq 0 ]; then
    echo "No tool name specified!"
    echo "Run build.sh -h for help"
    exit 1
fi

if [ $1 == "-h" ]; then
    echo "Utility for building docker containers from tools/"
    echo "Usage: build.sh [options]"
    echo 
    echo "Options:"
    echo "  all             build all tools including main image"
    echo "  <tool name>     builds specific tool"
    echo "                  Example: sh build.sh ericscript"
    exit 0
fi

if [ $1 == "all" ]; then
    for TOOL in */; do
        create_container `pwd`/$TOOL ${TOOL%?}
    done
    # Build main cotainer
    VERSION="$(cat ../nextflow.config | grep "container" | cut -d":" -f2 | cut -d "'" -f1)"
    CONTAINER_NAME=$PREFIX:$VERSION
    echo "Building [$CONTAINER_NAME]"
    docker build ../. -t $CONTAINER_NAME
else
    TOOL=$1
    TOOL_PATH="$(pwd)/$TOOL"
    if [ ! -d $TOOL_PATH ]; then
        echo "The tool doesn't exist"
        exit 1
    else
        create_container $TOOL_PATH
    fi
fi
