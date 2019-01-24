# nfcore/rnafusion

## nfcore/rnafusion version 1.0 -

Version 1.0 marks the first production release of this pipeline under the nf-core flag. The pipeline includes
additional help scripts to download references for fusion tools and Singularity images.

* Fusion gene detection tools:
  * `STAR-Fusion v1.5.0`
  * `Fusioncatcher v1.00 (commit hash: 045a8af)`
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