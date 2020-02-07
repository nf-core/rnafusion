
# ![nf-core/rnafusion](docs/images/nf-core-rnafusion_logo.png)

**RNA sequencing analysis pipeline with curated list of tools for detecting and visualizing fusion genes.**

[![Build Status](https://travis-ci.com/nf-core/rnafusion.svg?branch=master)](https://travis-ci.com/nf-core/rnafusion)
[![GitHub Actions CI Status](https://github.com/nf-core/rnafusion/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/rnafusion/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/rnafusion/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/rnafusion/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)
[![DOI](https://zenodo.org/badge/151721952.svg)](https://zenodo.org/badge/latestdoi/151721952)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/rnafusion.svg)](https://hub.docker.com/r/nfcore/rnafusion)

## Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

| Tool            | Single-end reads | CPU (recommended) | RAM (recommended) |
| --------------- |:----------------:|:-----------------:|:-----------------:|
| [Arriba](https://github.com/suhrig/arriba) | **No**           | >=16 cores        | ~30GB             |
| [EricScript](https://sites.google.com/site/bioericscript/getting-started) | **No**           | >=16 cores        | ~30GB             |
| [FusionCatcher](https://github.com/ndaniel/fusioncatcher) | Yes              | >=16 cores        | ~64GB             |
| [fusion-report](https://github.com/matq007/fusion-report) | -              | -        | -             |
| [Pizzly](https://github.com/pmelsted/pizzly) | **No**           | >=16 cores        | ~30GB             |
| [Squid](https://github.com/Kingsford-Group/squid) | **No**           | >=16 cores        | ~30GB             |
| [Star-Fusion](https://github.com/STAR-Fusion/STAR-Fusion) | Yes              | >=16 cores        | ~30GB             |
| [FusionInspector](https://github.com/FusionInspector/FusionInspector) | **No**           | >=16 cores        | ~30GB             |

For available parameters or help run:

```bash
nextflow run nf-core/rnafusion --help
```

## Quick Start

i. Install [`nextflow`](https://nf-co.re/usage/installation)

ii. Install one of [`docker`](https://docs.docker.com/engine/installation/), [`singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`conda`](https://conda.io/miniconda.html)

iii. Download the pipeline and test it on a minimal dataset with a single command

```bash
nextflow run nf-core/rnafusion --help
```

> Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile institute` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

iv. Start running your own analysis!

```bash
nextflow run nf-core/rnafusion -profile <profile> -c './example/custom-docker.config' --reads '*_R{1,2}.fastq.gz' --arriba --star_fusion --fusioncatcher --ericscript --pizzly --squid --arriba_vis --fusion_inspector
```

See [usage docs](docs/usage.md) for all of the available options when running the pipeline.

## Documentation

The nf-core/rnafusion pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
    * [Download references](docs/references.md)
    * [Local installation](https://nf-co.re/usage/local_installation)
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)

Use predefined configuration for desired Institution cluster provided at [nfcore/config](https://github.com/nf-core/configs) repository.

## Credits

This pipeline was originally written by Martin Proks ([@matq007](https://github.com/matq007)) in collaboration with Karolinska Institutet, SciLifeLab and University of Southern Denmark as a master thesis. This is a follow-up development started by Rickard Hammarén ([@Hammarn](https://github.com/Hammarn)). Special thanks goes to all supervisors: Teresita Díaz de Ståhl, PhD., Assoc. Prof.; Monica Nistér, MD, PhD; Maxime U Garcia PhD ([@MaxUlysse](https://github.com/MaxUlysse)); Szilveszter Juhos ([@szilvajuhos](https://github.com/szilvajuhos)); Phil Ewels PhD ([@ewels](https://github.com/ewels)) and Lars Grøntved, PhD., Assoc. Prof.

## Tool References

* **STAR-Fusion: Fast and Accurate Fusion Transcript Detection from RNA-Seq**
Brian Haas, Alexander Dobin, Nicolas Stransky, Bo Li, Xiao Yang, Timothy Tickle, Asma Bankapur, Carrie Ganote, Thomas Doak, Natalie Pochet, Jing Sun, Catherine Wu, Thomas Gingeras, Aviv Regev
bioRxiv 120295; doi: [https://doi.org/10.1101/120295](https://doi.org/10.1101/120295)
* D. Nicorici, M. Satalan, H. Edgren, S. Kangaspeska, A. Murumagi, O. Kallioniemi, S. Virtanen, O. Kilkku, **FusionCatcher – a tool for finding somatic fusion genes in paired-end RNA-sequencing data**, bioRxiv, Nov. 2014,
[DOI:10.1101/011650](http://dx.doi.org/10.1101/011650)
* Benelli M, Pescucci C, Marseglia G, Severgnini M, Torricelli F, Magi A. **Discovering chimeric transcripts in paired-end RNA-seq data by using EricScript**. Bioinformatics. 2012; 28(24): 3232-3239.
* **Fusion detection and quantification by pseudoalignment**
Páll Melsted, Shannon Hateley, Isaac Charles Joseph, Harold Pimentel, Nicolas L Bray, Lior Pachter, bioRxiv 166322; doi: [https://doi.org/10.1101/166322](https://doi.org/10.1101/166322)
* **SQUID: transcriptomic structural variation detection from RNA-seq** Cong Ma, Mingfu Shao and Carl Kingsford, Genome Biology, 2018, doi: [https://doi.org/10.1186/s13059-018-1421-5](https://doi.org/10.1186/s13059-018-1421-5)
* **Fusion-Inspector** download: [https://github.com/FusionInspector](https://github.com/FusionInspector)
* **fusion-report** download: [https://github.com/matq007/fusion-report](https://github.com/matq007/fusion-report); doi: [https://doi.org/10.5281/zenodo.3520171](https://doi.org/10.5281/zenodo.3520171)
* **FastQC** download: [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* **MultiQC** Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. [https://doi.org/10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354) Download: [https://multiqc.info/](https://multiqc.info/)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on [Slack](https://nfcore.slack.com/channels/rnafusion) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citation

If you use  nf-core/rnafusion for your analysis, please cite it using the following doi: [10.5281/zenodo.151721952](https://zenodo.org/badge/latestdoi/151721952)

You can cite the `nf-core` pre-print as follows:  

> Ewels PA, Peltzer A, Fillinger S, Alneberg JA, Patel H, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. **nf-core: Community curated bioinformatics pipelines**. *bioRxiv*. 2019. p. 610741. [doi: 10.1101/610741](https://www.biorxiv.org/content/10.1101/610741v1).

[![Barntumörbanken](docs/images/BTB_logo.png)](https://ki.se/forskning/barntumorbanken-0) | [![SciLifeLab](docs/images/SciLifeLab_logo.png)](https://scilifelab.se)
:-:|:-:
[![National Genomics Infrastructure](docs/images/NGI_logo.png)](https://ngisweden.scilifelab.se/) | [![University of Southern Denmark](docs/images/SDU_logo.png)](https://www.sdu.dk/da)
