#!/bin/bash
singularity pull --name rnafusion_star-fusion_v1.0dev.img docker://nfcore/rnafusion:star-fusion_v1.0dev
singularity pull --name rnafusion_fusioncatcher_v1.0dev.img docker://nfcore/rnafusion:fusioncatcher_v1.0dev
singularity pull --name rnafusion_fusion-inspector_v1.0dev.img docker://nfcore/rnafusion:fusion-inspector_v1.0dev
singularity pull --name rnafusion_ericscript_v1.0dev.img docker://nfcore/rnafusion:ericscript_v1.0dev
singularity pull --name rnafusion_pizzly_v1.0dev.img docker://nfcore/rnafusion:pizzly_v1.0dev
singularity pull --name rnafusion_squid_v1.0dev.img docker://nfcore/rnafusion:squid_v1.0dev
singularity pull --name rnafusion_summary-report_v1.0dev.img docker://nfcore/rnafusion:summary-report_v1.0dev
singularity pull --name rnafusion_dev.img docker://nfcore/rnafusion:dev
