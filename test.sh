#!/bin/bash

# Test tools only
nextflow run nf-core/rnafusion -profile munin -c /data1/references/rnafusion/1.1.0/munin-singularity.config --reads "/data1/references/rnafusion/1.1.0/reads_{1,2}.fq.gz" --arriba --star_fusion --fusioncatcher --ericscript --pizzly --squid --debug

# Test batch mode
nextflow run nf-core/rnafusion -profile munin -c /data1/references/rnafusion/1.1.0/munin-singularity.config --reads "/data1/references/rnafusion/1.1.0/batch_data/new_reads/*_R{1,2}.fq.gz" --arriba --star_fusion --fusioncatcher --ericscript --pizzly --squid

# Test the whole pipeline
nextflow run nf-core/rnafusion -profile munin -c /data1/references/rnafusion/1.1.0/munin-singularity.config --reads "/data1/references/rnafusion/1.1.0/reads_{1,2}.fq.gz" --arriba --star_fusion --fusioncatcher --ericscript --pizzly --squid --arriba_vis --fusion_inspector

# Test the whole pipeline in batch mode
nextflow run nf-core/rnafusion -profile munin -c /data1/references/rnafusion/1.1.0/munin-singularity.config --reads "/data1/references/rnafusion/1.1.0/batch_data/new_reads/*_R{1,2}.fq.gz" --arriba --star_fusion --fusioncatcher --ericscript --pizzly --squid --arriba_vis --fusion_inspector