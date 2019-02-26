#!/bin/bash

if [ $# -eq 0 ]; then
    echo "No output path specified"
    echo "Run download-singularity-img.sh -h for help"
    exit 1
fi

if [ $1 == "-h" ]; then
    echo "Utility for downloading singularity images from DockerHub"
    echo
    echo "Usage: download-singularity-img.sh [output directory]"
    echo "Example: sh download-singularity-img.sh /path"
    echo
    exit 0
elif [ -d "$1" ]; then
    cwd=$1/rnafusion_containers_v$RNAFUSION
    echo "Creating folder $cwd"
    mkdir $cwd && cd $cwd
    echo "Pulling images ..."
    
    RNAFUSION="$(cat ../nextflow.config | grep "container" | cut -d":" -f2 | cut -d "'" -f1)"
    STAR_FUSION="$(cat ../nextflow.config | grep "star_fusion_version" | cut -d"=" -f2 | tr -d \' | tr -d ' ')"
    FUSIONCATCHER="$(cat ../nextflow.config | grep "fusioncatcher_version" | cut -d"=" -f2 | tr -d \' | tr -d ' ')"
    ERICSCRIPT="$(cat ../nextflow.config | grep "ericscript_version" | cut -d"=" -f2 | tr -d \' | tr -d ' ')" 
    PIZZLY="$(cat ../nextflow.config | grep "pizzly_version" | cut -d"=" -f2 | tr -d \' | tr -d ' ')"
    SQUID="$(cat ../nextflow.config | grep "squid_version" | cut -d"=" -f2 | tr -d \' | tr -d ' ')" 
    FUSION_INSPECTOR="$(cat ../nextflow.config | grep "fusion_inspector_version" | cut -d"=" -f2 | tr -d \' | tr -d ' ')"
    
    singularity pull --name "rnafusion_v${RNAFUSION}.img" docker://nfcore/rnafusion:${RNAFUSION}
    singularity pull --name "rnafusion_star-fusion_v${STAR_FUSION}.img" docker://nfcore/rnafusion:star-fusion_v${STAR_FUSION}
    singularity pull --name "rnafusion_fusioncatcher_v${FUSIONCATCHER}.img" docker://nfcore/rnafusion:fusioncatcher_v${FUSIONCATCHER}
    singularity pull --name "rnafusion_ericscript_v${ERICSCRIPT}.img" docker://nfcore/rnafusion:ericscript_v${ERICSCRIPT}
    singularity pull --name "rnafusion_pizzly_v${PIZZLY}.img" docker://nfcore/rnafusion:pizzly_v${PIZZLY}
    singularity pull --name "rnafusion_squid_v${SQUID}.img" docker://nfcore/rnafusion:squid_v${SQUID}
    singularity pull --name "rnafusion_fusion-inspector_v${FUSION_INSPECTOR}.img" docker://nfcore/rnafusion:fusion-inspector_v${FUSION_INSPECTOR}
else
    echo "Path doesn't exists"
    exit 1
fi
