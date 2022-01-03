# nf-core/rnafusion: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.3.0dev nfcore/rnafusion - 2020/07/15

* Using official STAR-Fusion container [#160](https://github.com/nf-core/rnafusion/issues/160)

### Added

* Added social preview image [#107](https://github.com/nf-core/rnafusion/issues/107)
* Added support for GRCh37 genome assembly [#77](https://github.com/nf-core/rnafusion/issues/77)

### Changed

* Upgrade `fusion-report v2.1.2` to `fusion-report v2.1.3`
* Upgrade `fusion-report v2.1.1` to `fusion-report v2.1.2`
* Upgrade `fusion-report v2.1.0` to `fusion-report v2.1.1`
* Upgrade `Arriba v1.1.0` to `Arriba v1.2.0`
* Upgrade `fusion-report v2.0.2` to `fusion-report v2.1.0`

### Fixed

* Missing `strip-components` in `download-references.nf/star-fusion` [#148](https://github.com/nf-core/rnafusion/issues/148)
* Missing version prefix for cdna [#143](https://github.com/nf-core/rnafusion/issues/143)
* `samtools` missing header in empty file for FusionInspector [ref](https://github.com/STAR-Fusion/STAR-Fusion/issues/191)
* Removed `profile` from helper scripts [#139](https://github.com/nf-core/rnafusion/issues/139)
* Wrong url path for `Pfam-A.hmm.gz` [#140](https://github.com/nf-core/rnafusion/issues/140)

### Removed

* Removed `scripts/download-singularity-img.sh` and `download-singularity-img.nf` as they are not necessary any more

---

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

Initial release of nf-core/rnafusion, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
