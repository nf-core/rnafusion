# nfcore/rnafusion: Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/) and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## dev nfcore/rnafusion

### Added

* Added social preview image [#107](https://github.com/nf-core/rnafusion/issues/107)

## [1.1.0] nfcore/rnafusion - 2020/02/10

* Fusion gene detection tools:
  * `Arriba v1.1.0`
  * `Ericscript v0.5.5`
  * `Fusioncatcher v1.20`
  * `Pizzly v0.37.3`
  * `Squid v1.5`
  * `STAR-Fusion v1.6.0`
* Visualization tools:
  * `Arriba v1.1.0`
  * `FusionInspector v1.3.1`
* Other tools:
  * `fusion-report v2.0.1`
  * `FastQ v0.11.8`
  * `MultiQC v1.7`
  * `STAR aligner v2.7.0f`

### Added

* Added `Arriba 1.1.0` [#63](https://github.com/nf-core/rnafusion/issues/63)
* Added Batch mode [#54](https://github.com/nf-core/rnafusion/issues/54)

### Changed

* Updated examples and configurations
* Upgraded `fusion-report v1.0.0` to `fusion-report v2.0.1`
* Divided `running_tools` into fusion and visualization tools
* Updated `STAR` in `Squid`, `Fusion-Inspector` version to `2.7.0f`
* Upgraded `STAR-Fusion v1.5.0` to `STAR-Fusion v1.6.0` [#83](https://github.com/nf-core/rnafusion/issues/83)
* Parameter `igenomesIgnore` renamed to `igenome` [#81](https://github.com/nf-core/rnafusion/issues/81)
* Finished STAR-Fusion file renaming [#18](https://github.com/nf-core/rnafusion/issues/18)
* Updated logos
* Updated to nf-core `1.8` TEMPLATE

### Fixed

* iGenomes optional, but not really [#91](https://github.com/nf-core/rnafusion/issues/91)
* Updated `fusioncatcher` to latest `1.20` version also solving [#95](https://github.com/nf-core/rnafusion/issues/95)

### Removed

* Variables `pizzly_fasta` and `pizzly_gtf` have been removed and replaced with `transcript` and `gtf`
* `Jenkisfile`, test configuration, pylintrc configuration
* Removed `igenomes.config` because the pipeline only supports `Ensembl` version

---

## [1.0.2] nfcore/rnafusion - 2019/05/13

### Changed

* Bumped nf-core template to 1.6 [#69](https://github.com/nf-core/rnafusion/pull/69)

### Fixed

* Fixed COSMIC parameters not wrapped in quotes [#75](https://github.com/nf-core/rnafusion/issues/75)
* Implemented output output for fusion tools [#72](https://github.com/nf-core/rnafusion/issues/72)
* Fixed reference download link for STAR-Fusion [#71](https://github.com/nf-core/rnafusion/issues/71)

---

## [1.0.1] nfcore/rnafusion - 2019/04/06

### Added

* Added support for extra parameters for tools STAR-Fusion, FusionCatcher and fusion-report
* Added example configuration for `singularity` and `docker`
* Added [fusion-report](https://github.com/matq007/fusion-report) into the stack [#62](https://github.com/nf-core/rnafusion/issues/62), [#55](https://github.com/nf-core/rnafusion/issues/55), [#53](https://github.com/nf-core/rnafusion/issues/53), [#51](https://github.com/nf-core/rnafusion/issues/51)
* Added nextflow helper script `download-singularity-img.nf`
* Added nextflow helper script `download-references.nf`
* Added `Jenkinsfile` for in-house testing

### Changed

* Updated installation of `FusionCatcher` (available now on bioconda)

### Fixed

* Fixed empty symlinks (`input.X`) in fusion-report [#68](https://github.com/nf-core/rnafusion/issues/68)
* Fixed FASTA issues [#60](https://github.com/nf-core/rnafusion/issues/60)
* Fixed centralized nf-core/config [#64](https://github.com/nf-core/rnafusion/issues/64)
* Fixed `scrape_software_versions.py` to parse tools versions correctly [#65](https://github.com/nf-core/rnafusion/issues/65)

### Removed

* Removed `Singularity`

---

## [1.0] nfcore/rnafusion - 2018/02/14

Version 1.0 marks the first production release of this pipeline under the nf-core flag.
The pipeline includes additional help scripts to download references for fusion tools and Singularity images.

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

---

## [0.1] SciLifeLab/NGI-RNAfusion (ARCHIVED) - 2018/10/05

Initial release of NGI-RNAfusion, created with the [nf-core](http://nf-co.re/) template.
Source code can be found at [SciLifeLab/NGI-RNAfusion](https://github.com/SciLifeLab/NGI-RNAfusion).
The solution works with Docker and Singularity.

* Tools:
  * STAR-Fusion
  * Fusioncatcher
  * FusionInspector
  * Custom tool for fusion comparison - generates intersection of detected fusion genes from all tools
