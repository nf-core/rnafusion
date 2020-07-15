# ![nf-core/rnafusion](docs/images/nf-core-rnafusion_logo.png)

**RNA sequencing analysis pipeline with curated list of tools for detecting and visualizing fusion genes.**

[![GitHub Actions CI Status](https://github.com/nf-core/rnafusion/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/rnafusion/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/rnafusion/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/rnafusion/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)
[![DOI](https://zenodo.org/badge/151721952.svg)](https://zenodo.org/badge/latestdoi/151721952)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/rnafusion.svg)](https://hub.docker.com/r/nfcore/rnafusion)

## Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

> The pipeline **requires** >=16 CPU cores and >=30GB RAM

| Tool                                                                      |  Single-end reads  |  Version |
| ------------------------------------------------------------------------- | :----------------: | :------: |
| [Arriba](https://github.com/suhrig/arriba)                                |         :x:        |  `1.2.0` |
| [EricScript](https://sites.google.com/site/bioericscript/getting-started) |         :x:        |  `0.5.5` |
| [FusionCatcher](https://github.com/ndaniel/fusioncatcher)                 | :white_check_mark: |  `1.20`  |
| [Fusion-Inspector](https://github.com/FusionInspector/FusionInspector)    |         :x:        |  `2.2.1` |
| [fusion-report](https://github.com/matq007/fusion-report)                 |          -         |  `2.1.3` |
| [Pizzly](https://github.com/pmelsted/pizzly)                              |         :x:        | `0.37.3` |
| [Squid](https://github.com/Kingsford-Group/squid)                         |         :x:        |   `1.5`  |
| [Star-Fusion](https://github.com/STAR-Fusion/STAR-Fusion)                 | :white_check_mark: |  `1.8.1` |

For available parameters or help run:

```bash
nextflow run nf-core/rnafusion --help
```

## Quick Start

i. Install [`nextflow`](https://nf-co.re/usage/installation)

ii. Install either [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) for full pipeline reproducibility (please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))

iii. Download references for all tools

```bash
nextflow run nf-core/rnafusion/download-references.nf -profile <docker/singularity/institute> \
  --download_all \
  --outdir <PATH> \
  --cosmic_usr <COSMIC_USER> --cosmic_passwd <COSMIC_PASSWD>
```

> Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

iv. Start running your own analysis!

```bash
nextflow run nf-core/rnafusion -profile <docker/singularity/institute> \
  --reads '*_R{1,2}.fastq.gz' \
  --genomes_base 'reference_path_from_above'
  --arriba --star_fusion --fusioncatcher --ericscript --pizzly --squid \
  --arriba_vis --fusion_inspector
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

This pipeline was originally written by Martin Proks ([@matq007](https://github.com/matq007)) in collaboration with Karolinska Institutet, SciLifeLab and University of Southern Denmark as a master thesis. This is a follow-up development started by Rickard Hammarén ([@Hammarn](https://github.com/Hammarn)).

Special thanks goes to all supervisors:

* [Assoc. Prof. Teresita Díaz de Ståhl, PhD](https://ki.se/en/onkpat/teresita-diaz-de-stahls-group)
* [MD. Monica Nistér, PhD](https://ki.se/en/onkpat/research-team-monica-nister)
* [Maxime U Garcia, PhD](https://github.com/MaxUlysse)
* [Szilveszter Juhos](https://github.com/szilvajuhos)
* [Phil Ewels, PhD](https://github.com/ewels)
* [Assoc. Prof. Lars Grøntved, PhD](https://portal.findresearcher.sdu.dk/en/persons/larsgr)

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

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).  
> ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)

[![Barntumörbanken](docs/images/BTB_logo.png)](https://ki.se/forskning/barntumorbanken-0) | [![SciLifeLab](docs/images/SciLifeLab_logo.png)](https://scilifelab.se)
:-:|:-:
[![National Genomics Infrastructure](docs/images/NGI_logo.png)](https://ngisweden.scilifelab.se/) | [![University of Southern Denmark](docs/images/SDU_logo.png)](https://www.sdu.dk/da)
