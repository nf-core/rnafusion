#!/bin/bash

get_tool_version() {
    echo $(cat nextflow.config | grep -m1 "$1" | cut -d"=" -f2 | tr -d \' | tr -d ' ')
}

if [ $# -eq 0 ]; then
    echo "No output path specified"
    echo "Run download-singularity-img.sh -h for help"
    exit 1
fi

if [ $1 == "-h" ]; then
    echo "Utility for downloading singularity images from DockerHub"
    echo
    echo "Usage: download-singularity-img.sh [output directory]"
    echo "Example: sh scripts/download-singularity-img.sh /path"
    echo
    exit 0
elif [ -d "$1" ]; then
    PREFIX="nf-core-rnafusion"
    ARRIBA=$(get_tool_version "arriba")
    ERICSCRIPT=$(get_tool_version "ericscript")
    FUSIONCATCHER=$(get_tool_version "fusioncatcher")
    FUSION_INSPECTOR=$(get_tool_version "fusion_inspector")
    PIZZLY=$(get_tool_version "pizzly")
    RNAFUSION="$(cat nextflow.config | grep -m1 "container" | cut -d":" -f2 | cut -d "'" -f1)"
    SQUID=$(get_tool_version "squid")
    STAR_FUSION=$(get_tool_version "star_fusion")

    cd $1 && echo "Pulling images ..."

    singularity pull --name "${PREFIX}-arriba_${ARRIBA}.img" docker://nfcore/rnafusion:arriba_${ARRIBA}
    singularity pull --name "${PREFIX}-ericscript_${ERICSCRIPT}.img" docker://nfcore/rnafusion:ericscript_${ERICSCRIPT}
    singularity pull --name "${PREFIX}-fusioncatcher_${FUSIONCATCHER}.img" docker://nfcore/rnafusion:fusioncatcher_${FUSIONCATCHER}
    singularity pull --name "${PREFIX}-fusion-inspector_${FUSION_INSPECTOR}.img" docker://nfcore/rnafusion:fusion-inspector_${FUSION_INSPECTOR}
    singularity pull --name "${PREFIX}-pizzly_${PIZZLY}.img" docker://nfcore/rnafusion:pizzly_${PIZZLY}
    singularity pull --name "${PREFIX}-${RNAFUSION}.img" docker://nfcore/rnafusion:${RNAFUSION}
    singularity pull --name "${PREFIX}-squid_${SQUID}.img" docker://nfcore/rnafusion:squid_${SQUID}
    singularity pull --name "${PREFIX}-star-fusion_${STAR_FUSION}.img" docker://nfcore/rnafusion:star-fusion_${STAR_FUSION}

else
    echo "Path doesn't exists"
    exit 1
fi