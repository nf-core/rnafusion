# ![nf-core/rnafusion](https://raw.githubusercontent.com/nf-core/rnafusion/master/docs/images/rnafusion_logo.png)

[![Build Status](https://travis-ci.org/nf-core/rnafusion.svg?branch=master)](https://travis-ci.org/nf-core/rnafusion)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/rnafusion.svg)](https://hub.docker.com/r/nfcore/rnafusion)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)

### Introduction

**nfcore/rnafusion** uses RNA-seq data to detect fusions genes.

The workflow processes raw single-read (limited tools available) or paired-end data from FastQ input ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)), detect fusion genes ([STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion), [Fusioncatcher](https://github.com/ndaniel/fusioncatcher), [Ericscript](https://sites.google.com/site/bioericscript/), [Pizzly](https://github.com/pmelsted/pizzly), [Squid](https://github.com/Kingsford-Group/squid)), visualizes the fusions ([FusionInspector](https://github.com/FusionInspector/FusionInspector)), performs quality-control on the results ([MultiQC](http://multiqc.info)) and finally generates custom summary report.

**Note: Make sure to read the installation guide before running the pipeline.**

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.

### Documentation
The nf-core/rnafusion pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Tools](docs/tools.md)
    * [Local installation](docs/configuration/local.md)
    * [Swedish UPPMAX cluster](docs/configuration/uppmax.md)
    * [Adding your own system](docs/configuration/adding_your_own.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

### Final summary output

![Final summary report](docs/images/example-summary-report.png)