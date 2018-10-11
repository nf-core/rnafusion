#!/bin/bash

# Run pipeline
nextflow run NGI-RNAfusion -profile test --reads "/test-data/test_{1,2}.fastq.gz" --genome R64-1-1 --star_fusion

# Test per tool
# cp -r /results/* NGI-RNAfusion/results
