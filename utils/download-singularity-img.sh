#!/bin/bash
VERSION="$(cat ../nextflow.config | grep "container_version" | cut -d"=" -f2 | tr -d "\'" | tr -d " ")"

singularity pull --name rnafusion_star-fusion_v${VERSION}.img docker://nfcore/rnafusion:star-fusion_v${VERSION}
singularity pull --name rnafusion_fusioncatcher_v${VERSION}.img docker://nfcore/rnafusion:fusioncatcher_v${VERSION}
singularity pull --name rnafusion_fusion-inspector_v${VERSION}.img docker://nfcore/rnafusion:fusion-inspector_v${VERSION}
singularity pull --name rnafusion_ericscript_v${VERSION}.img docker://nfcore/rnafusion:ericscript_v${VERSION}
singularity pull --name rnafusion_pizzly_v${VERSION}.img docker://nfcore/rnafusion:pizzly_v${VERSION}
singularity pull --name rnafusion_squid_v${VERSION}.img docker://nfcore/rnafusion:squid_v${VERSION}
singularity pull --name rnafusion_summary-report_v${VERSION}.img docker://nfcore/rnafusion:summary-report_v${VERSION}
singularity pull --name rnafusion-${VERSION}.img docker://nfcore/rnafusion:${VERSION}
