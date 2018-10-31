#!/bin/bash

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
        VERSION=$(grep $(head -1 $TOOL_PATH/environment.yml | cut -d' ' -f2)'=' $TOOL_PATH/environment.yml | cut -d'=' -f2)
        CONTAINER_NAME=nfcore/rnafusion:${TOOL}_v${VERSION}
        docker build $TOOL_PATH -t $TOOL
        docker tag $TOOL $CONTAINER_NAME
        docker push $CONTAINER_NAME
    fi
fi
