#!/bin/bash

VERSION="$(cat ../nextflow.config | grep "container_version" | cut -d"=" -f2 | tr -d "\'" | tr -d " ")"
if [ $# -eq 0 ]; then
    echo "No tool name specified!"
    exit 1
else
    TOOL=$1
    TOOL_PATH="$(pwd)/$TOOL"
    if [ ! -d $TOOL_PATH ]; then
        echo "The tool doesn't exist"
        exit 1
    else
        CONTAINER_NAME=nfcore/rnafusion:${TOOL}_v${VERSION}
        docker build $TOOL_PATH -t $CONTAINER_NAME
    fi
fi
