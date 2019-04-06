# nfcore/rnafusion

## nfcore/rnafusion version 1.0.1 -

* Added example configuration for `singularity` and `docker`
* Removed `Singularity`
* Fixed FASTA issues [#60](https://github.com/nf-core/rnafusion/issues/60)
* Fixed centralized nf-core/config [#64](https://github.com/nf-core/rnafusion/issues/64)
* Added [fusion-report](https://github.com/matq007/fusion-report) into the stack [#62](https://github.com/nf-core/rnafusion/issues/62)[#55](https://github.com/nf-core/rnafusion/issues/55)[#53](https://github.com/nf-core/rnafusion/issues/53)[#51](https://github.com/nf-core/rnafusion/issues/51)
* Added nextflow helper script `download-singularity-img.nf`
* Added nextflow helper script `download-references.nf`
* Fixed `scrape_software_versions.py` to parse tools versions correctly
* Added `Jenkinsfile` for in-house testing
* Updated installation of `FusionCatcher` (available now on bioconda)

## nfcore/rnafusion version 1.0 - 2018/02/14

Version 1.0 marks the first production release of this pipeline under the nf-core flag. The pipeline includes
additional help scripts to download references for fusion tools and Singularity images.

* Fusion gene detection tools:
  * `STAR-Fusion v1.5.0`
  * `Fusioncatcher v1.00`
  * `Ericscript v0.5.5`
  * `Pizzly v0.37.3`
  * `Squid v1.5`
* Visualization tools:
  * `FusionInspector v1.3.1`
* Other tools:
  * `Summary report`
  * `FastQ v0.11.8`
  * `MultiQC v1.7`
  * `FusionGDB updated 2019/01/23`

## SciLifeLab/NGI-RNAfusion version 0.1 (ARCHIVED) - 2018/10/05

Initial release of NGI-RNAfusion, created with the [nf-core](http://nf-co.re/) template. Source code can be found
at [SciLifeLab/NGI-RNAfusion](https://github.com/SciLifeLab/NGI-RNAfusion). The solution works with Docker and Singularity.

* Tools:
  * STAR-Fusion
  * Fusioncatcher
  * FusionInspector
  * Custom tool for fusion comparison - generates intersection of detected fusion genes from all tools