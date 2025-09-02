# nf-core/rnafusion: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v4.0.0dev - [date]

### Added

- Normalized gene expression calculated [#488](https://github.com/nf-core/rnafusion/pull/488)
- Primary assembly now used as main reference genome FASTA file, as recommended by the STAR manual [#488](https://github.com/nf-core/rnafusion/pull/488)
- Use of only ensembl GTF file, not chr.gtf file as GTF reference file [#488](https://github.com/nf-core/rnafusion/pull/488)
- Add nf-test to local module: `ENSEMBL_DOWNLOAD` [#539](https://github.com/nf-core/rnafusion/pull/539)
- Add nf-test to local module: `HGNC_DOWNLOAD` [#540](https://github.com/nf-core/rnafusion/pull/540)
- Add nf-test to local subworkflow: `STRINGTIE_WORKFLOW` [#541](https://github.com/nf-core/rnafusion/pull/541)
- Option to avoid using COSMIC (for example in the case of clinical use) [#547](https://github.com/nf-core/rnafusion/pull/547)
- Add nf-test to nf-core module: `PICARD_COLLECTRNASEQMETRICS` and update module [#551](https://github.com/nf-core/rnafusion/pull/551)
- Add `--skip_vcf` boolean parameter to skip vcf file generation [#554](https://github.com/nf-core/rnafusion/pull/554)
- Add nf-test to local module: `FUSIONREPORT_DOWNLOAD` [#560](https://github.com/nf-core/rnafusion/pull/560)
- Add nf-test to local subworkflow: `QC_WORKFLOW` [#568](https://github.com/nf-core/rnafusion/pull/568)
- Add nf-test to local subworkflow: `TRIM_WORKFLOW` [#572](https://github.com/nf-core/rnafusion/pull/572)
- Add nf-test to local module: `FUSIONREPORT_DETECT`. Improve `FUSIONREPORT_DOWNLOAD` module [#577](https://github.com/nf-core/rnafusion/pull/577)
- Add nf-test to local subworkflow: `ARRIBA_WORKFLOW` [#578](https://github.com/nf-core/rnafusion/pull/578)
- Add nf-test to local module: `STARFUSION_BUILD`. [#585](https://github.com/nf-core/rnafusion/pull/585)
- Add nf-test to local module: `STARFUSION_DETECT`. [#586](https://github.com/nf-core/rnafusion/pull/586)
- Added a new module `CTATSPLICING_STARTOCANCERINTRONS` and a new parameter `--ctatsplicing`. This options creates reports on cancer splicing aberrations and requires one or both of `--arriba` and `--starfusion` to be given. [#587](https://github.com/nf-core/rnafusion/pull/587)
- Add parameter `--references_only` when no data should be analysed, but only the references should be built [#505](https://github.com/nf-core/rnafusion/pull/505)
- Add nf-test to local subworkflow: `FUSIONCATCHER_WORKFLOW` [#591](https://github.com/nf-core/rnafusion/pull/591)
- Add nf-test to local subworkflow: `STARFUSION_WORKFLOW`. [#597](https://github.com/nf-core/rnafusion/pull/597)
- Add nf-test to local module: `FUSIONINSPECTOR`. [#601](https://github.com/nf-core/rnafusion/pull/601)
- Added `CTATSPLICING_PREPGENOMELIB` to update the starfusion genome library directory with a cancer splicing index. [#610](https://github.com/nf-core/rnafusion/pull/610)
- Add nf-test to local subworkflow: `FUSIONREPORT_WORKFLOW`. [#607](https://github.com/nf-core/rnafusion/pull/607)
- Add nf-test to local module: `ARRIBA_VISUALISATION`. [#625](https://github.com/nf-core/rnafusion/pull/625)
- Added the following fields to the samplesheet [#647](https://github.com/nf-core/rnafusion/pull/647):
  - `bam`: A BAM file aligned with STAR, it's the responsibility of the pipeline user to make sure this file has been correctly called.
  - `bai`: The index of the BAM file, this is not required when a `bam` file has been given but can increase the pipeline speed a bit.
  - `cram`: A CRAM file aligned with STAR, it's the responsibility of the pipeline user to make sure this file has been correctly called.
  - `crai`: The index of the CRAM file, this is not required when a `cram` file has been given but can increase the pipeline speed a bit.
  - `junctions`: A file containing the junctions determined by STAR (needed by `starfusion` and `ctatsplicing`)
  - `splice_junctions` A file containing the splice junctions determined by STAR (needed by `ctatsplicing`)
- Added `--fusioncatcher_download_link`. [#650](https://github.com/nf-core/rnafusion/pull/650)
- Added the `seq_platform` and `seq_center` fields to the samplesheet. These values can be used to overwrite the value of `--seq_platform` and `--seq_center` on a sample-by-sample basis [#654](https://github.com/nf-core/rnafusion/pull/654)
- Moved strandedness determination for picard collectrnaseqmetrics into modules.config [#658](https://github.com/nf-core/rnafusion/pull/658)
- Using option `--human_gencode_filter` while building STARFusion references for human species [#657](https://github.com/nf-core/rnafusion/pull/657)
- Added `--star_limit_bam_sort_ram` to set the maximum memory for sorting the BAM files in STAR. [#668](https://github.com/nf-core/rnafusion/pull/668)
- Update STAR-Fusion to 1.15.0 and update nf-core modules [#664](https://github.com/nf-core/rnafusion/pull/664)
- Update fusionreport 4.0.0 that replaces the score of a fusion with a Fusion Indication Index [#667](https://github.com/nf-core/rnafusion/pull/667)
- Added a extra trimming step for fusioncatcher in case a different read-length is wished for this tool only [#674](https://github.com/nf-core/rnafusion/pull/674)
- Added support for references hosted on the nf-core AWS S3 bucket. [#717](https://github.com/nf-core/rnafusion/pull/717)
- Added `--dfam_hmm`, `--dfam_h3f`, `--dfam_h3i`, `--dfam_h3m`, `--dfam_h3p`, `--pfam_file` and `--annot_filter_url` parameters to allow use of custom files in `STARFUSION_BUILD` module [#709](https://github.com/nf-core/rnafusion/pull/709)
- Add nf-test to local module: `VCF_COLLECT`. [#745](https://github.com/nf-core/rnafusion/pull/745)
- Added nf-test for local subworkflow: `FUSIONINSPECTOR_WORKFLOW`. [#753](https://github.com/nf-core/rnafusion/pull/753)

### Changed

- Updated modules and migrated non-specific modules to nf-core/modules [#484](https://github.com/nf-core/rnafusion/pull/484)
- Updated to nf-core/tools 3.0.2 [#504](https://github.com/nf-core/rnafusion/pull/504)
- Remove local module `RRNA_TRANSCRIPTS` (replaced by nf-core module) [#541](https://github.com/nf-core/rnafusion/pull/541)
- Allow fastq files without a dot before .fn(.gz)/.fastq(.gz) files [#548](https://github.com/nf-core/rnafusion/pull/548)
- Remove double nested folder introduced in [#577](https://github.com/nf-core/rnafusion/pull/577), [#581](https://github.com/nf-core/rnafusion/pull/581)
- Use docker.io and galaxy containers for fusioncatcher and starfusion (incl. fusioninspector) instead of wave as they are not functional on wave [#588](https://github.com/nf-core/rnafusion/pull/588)
- Update STAR-Fusion to 1.14 [#588](https://github.com/nf-core/rnafusion/pull/588)
- Use "-genePredExt -geneNameAsName2 -ignoreGroupsWithoutExons" (to mimic gms/tomte) for GTF_TO_REFFLAT [#505](https://github.com/nf-core/rnafusion/pull/505)
- Integrate reference building in the main workflow [#505](https://github.com/nf-core/rnafusion/pull/505)
- Move from ensembl to gencode base [#505](https://github.com/nf-core/rnafusion/pull/505)
- Update from ensembl 102 to gencode 46 default references [#505](https://github.com/nf-core/rnafusion/pull/505)
- Update`FUSIONINSPECTOR` to v2.10.0. [#601](https://github.com/nf-core/rnafusion/pull/601)
- Remove local module `STARFUSION_DOWNLOAD` [#598](https://github.com/nf-core/rnafusion/pull/598)
- Fix error message when parameter outdir is missing [#611](https://github.com/nf-core/rnafusion/pull/611)
- Updated documentation for fusion-report score calculation to reflect 80/20 weight distribution between thttps://github.com/nf-core/rnafusion/pull/633ool detection and database hits [#620](https://github.com/nf-core/rnafusion/pull/620)
- The STAR alignment now only runs once instead of multiple times when using `--arriba` and `--starfusion` [#633](https://github.com/nf-core/rnafusion/pull/633)
- The `--cram` parameter has been converted to a boolean value instead of the comma-separated list of values. Use this parameter if you also want to create CRAM files from the BAM files created with STAR [#633](https://github.com/nf-core/rnafusion/pull/633)
- Update htslib and samtools version in `star/align` to 1.21 [#634](https://github.com/nf-core/rnafusion/pull/634)
- Updated Zenodo DOI in badge [#639](https://github.com/nf-core/rnafusion/pull/639)
- `--run_fusioncatcher` back to `fusioncatcher` [#641](https://github.com/nf-core/rnafusion/pull/641)
- Removed `--fastp_trim`, `--arriba`, `--ctatsplicing`, `--fusioncatcher`, `--starfusion`, `--stringtie` and `--all` and replaced these parameters with the `--tools` parameter. This parameter takes a comma-delimited list of tool names to run for the pipeline and is a required parameter. Following tools are supported by this parameter [#645](https://github.com/nf-core/rnafusion/pull/645):
  - `arriba`
  - `ctatsplicing`
  - `fusioncatcher`
  - `starfusion`
  - `stringtie`
  - `fusionreport`
  - `fastp`
  - `salmon`
  - `fusioninspector`
  - `all` => This will automatically run all of the above tools
- Updated all `wget` containers to `conda-forge::wget=1.21.4` [#655](https://github.com/nf-core/rnafusion/pull/655)
- Using FusionInspector abridged output (that now contains coding effects) for downstream analysis [#669](https://github.com/nf-core/rnafusion/pull/669)
- Replaced local `FUSIONREPORT` and `FUSIONREPORT_DOWNLOAD` with `FUSIONREPORT_DETECT` and `FUSIONREPORT_DOWNLOAD` from nf-core [#703](https://github.com/nf-core/rnafusion/pull/703)
- Updated arriba to v2.5.0 [#693](https://github.com/nf-core/rnafusion/pull/693)
- Update logo [#715](https://github.com/nf-core/rnafusion/pull/715)
- Replaced local `STARFUSION_BUILD` for module from nf-core [#709](https://github.com/nf-core/rnafusion/pull/709)
- Modified `test_build` profile to use a reduced version of Pfam and Dfam files [#733](https://github.com/nf-core/rnafusion/pull/733)
- Changed local `ARRIBA_VISUALIZATION`, `CTATSPLICING_STARTOCANCERINTRONS`, `CTATSPLICING_PREPGENOMELIB`, `FUSIONINSPECTOR`, `STARFUSION_DETECT` for its nf-core module versions [#740](https://github.com/nf-core/rnafusion/pull/740)
- Replaced local subworkflow `TRIM_WORKFLOW` for its nf-core subworkflow equivalent `FASTQ_FASTQC_UMITOOLS_FASTP` [#752](https://github.com/nf-core/rnafusion/pull/752)
- Changed local `FASTQ_ALIGN_STAR` for subworkflow from nf-core [#756](https://github.com/nf-core/rnafusion/pull/756)

### Fixed

- Fixed some Nextflow run-commands in the docs [#491](https://github.com/nf-core/rnafusion/pull/491)
- Fixed bug when trying to build indices behind a proxy and wget was unable to download arriba indices [#495](https://github.com/nf-core/rnafusion/issues/495)
- Fixed bug in `FUSIONREPORT_DOWNLOAD` when building references with `--no_cosmic parameter` [#555](https://github.com/nf-core/rnafusion/issues/555)
- Refactor structure in `FUSIONREPORT_DOWNLOAD` to use cosmic credentials in `ext.args` [#556](https://github.com/nf-core/rnafusion/issues/556)
- Fixed bug in nf-core `RRNATRANSCRIPTS` module [#563](https://github.com/nf-core/rnafusion/issues/563)
- Fixed bug in `GFFREAD` that caused output `gffread_fasta` not being produced [#565](https://github.com/nf-core/rnafusion/issues/565)
- Fixed bug in `FUSIONCATCHER_DOWNLOAD` that caused an error when running with singularity profile [#573](https://github.com/nf-core/rnafusion/issues/573)
- Fixed missing script `gtf2bed` which caused local module `GET_RRNA_TRANSCRIPTS` to fail [#602](https://github.com/nf-core/rnafusion/issues/602)
- Fixed the codebase to be compatible with the Nextflow language server [#634](https://github.com/nf-core/rnafusion/pull/634)
- Updated the input validation to be more strict. This will prevent more errors down the line in the pipeline [#640](https://github.com/nf-core/rnafusion/pull/640)
- The `FUSIONINSPECTOR` process will no longer fail when no fusions have been found. [#651](https://github.com/nf-core/rnafusion/pull/651)
- Fixed STAR-Fusion build would fail when downloading Pfam and Dfam resources behind SSL. [#653](https://github.com/nf-core/rnafusion/pull/653)
- Fix bug in argument handling for FusionInspector [#669](https://github.com/nf-core/rnafusion/pull/669)
- Upgrade STAR-Fusion to 1.15.1 to solve problem building STAR-Fusion references with `--human_gencode_filter` [#683](https://github.com/nf-core/rnafusion/pull/683)
- Fix missing memory unit for fusioncatcher [#674](https://github.com/nf-core/rnafusion/pull/674)
- Fix fusioncatcher download link [#693](https://github.com/nf-core/rnafusion/pull/693)
- Fix fusionreport singularity container [#713](https://github.com/nf-core/rnafusion/pull/713)
- Fix CTAT-SPLICING output when no cancer introns were found [#722](https://github.com/nf-core/rnafusion/pull/722)
- Update VCF_COLLECT script to adapt to transcript_version not being an entry in fusioninspector gtf anymore [#726](https://github.com/nf-core/rnafusion/pull/726)
- Fix rRNA detection and make it more customizable with nf-core modules [#736](https://github.com/nf-core/rnafusion/pull/736)
- Fix rRNA detection in GTF using `transcript_type` [#749](https://github.com/nf-core/rnafusion/pull/749)
- Fixed nf-test for `QC_WORKFLOW` subworkflow [#756](https://github.com/nf-core/rnafusion/pull/756)

### Removed

- Remove fusionGDB from documentation and fusion-report download stubs [#503](https://github.com/nf-core/rnafusion/pull/503)
- Removed test-build as reference building gets integrated in the main workflow [#505](https://github.com/nf-core/rnafusion/pull/505)
- Removed parameter `--build_references`
- Removed fusioncatcher build [#650](https://github.com/nf-core/rnafusion/pull/650)
- Removed subworkflow with less than two modules: `ARRIBA_WORKFLOW` [#692](https://github.com/nf-core/rnafusion/pull/692)
- Removed subworkflow with less than two modules: `STARFUSION_WORKFLOW` [#707](https://github.com/nf-core/rnafusion/pull/707)
- Removed subworkflow with less than two modules: `CTATSPLICING_WORKFLOW` [#704](https://github.com/nf-core/rnafusion/pull/704)
- Removed subworkflow with less than two modules: `FUSIONREPORT_WORKFLOW` [#721](https://github.com/nf-core/rnafusion/pull/721)
- Removed local module `GET_RRNA_TRANSCRIPTS` [#736](https://github.com/nf-core/rnafusion/pull/736)
- Removed old unused parameters `params.download_refs` and `params.fusioncatcher_download_link` [#752](https://github.com/nf-core/rnafusion/pull/752)

### Parameters

| Old parameter         | New parameter                   |
| --------------------- | ------------------------------- |
|                       | `--no_cosmic`                   |
| `--build_references`  | `--references_only`             |
| `--fastp_trim`        | `--tools fastp`                 |
| `--arriba`            | `--tools arriba`                |
| `--run_fusioncatcher` | `--tools fusioncatcher`         |
| `--starfusion`        | `--tools starfusion`            |
| `--stringtie`         | `--tools stringtie`             |
| `--all`               | `--tools all`                   |
|                       | `--fusioncatcher_download_link` |
|                       | `--trim_tail_fusioncatcher`     |
|                       | `--save_trimmed_fail`           |
|                       | `--save_merged`                 |
|                       | `--min_trimmed_reads`           |
|                       | `--trim_tail_fusioncatcher`     |

## v3.0.2 - [2024-04-10]

### Added

### Changed

- Update to nf-tools 2.11.1 [#457](https://github.com/nf-core/rnafusion/pull/457)
- Update picard collectrnaseqmetrics memory requirements to 0.8x what is provided [#474](https://github.com/nf-core/rnafusion/pull/474)

### Fixed

- Fix bug when using parameter "whitelist" [#466](https://github.com/nf-core/rnafusion/pull/466)
- Fix VCF_COLLECT handling when a tool is absent from FUSIONREPORT report [#458](https://github.com/nf-core/rnafusion/pull/458)
- Fix VCF_COLLECT when fusioninspector output is empty but fusionreport is not [#465](https://github.com/nf-core/rnafusion/pull/465)
- Fix VCF_COLLECT bug [#481](https://github.com/nf-core/rnafusion/pull/481)
- Fix conda package for starfusion/detect[#482](https://github.com/nf-core/rnafusion/pull/482)
- Fix logical gate so when stringtie should run but not starfusion, starfusion will not run[#482](https://github.com/nf-core/rnafusion/pull/482)

### Removed

## v3.0.1 - [2023-11-29]

### Added

### Changed

- Python3 explicit in vcf_collect [#452](https://github.com/nf-core/rnafusion/pull/452)

### Fixed

- software-version.yml and in general version track-keeping was incomplete [#451](https://github.com/nf-core/rnafusion/pull/451)

### Removed

## v3.0.0 - [2023-11-27]

### Added

- Add picard CollectInsertSizeMetrics to QC workflow [#408](https://github.com/nf-core/rnafusion/pull/408)
- Build CRAM index in the same directory as CRAM files for Arriba and STAR-Fusion [#427](https://github.com/nf-core/rnafusion/pull/427)

### Changed

- Replace PICARD_MARKDUPLICATES with GATK4_MARKDUPLICATES [#409](https://github.com/nf-core/rnafusion/pull/409)
- Removed `--fusioninspector_filter` and `--fusionreport_filter` in favor of `--tools_cutoff` (default = 1, no filters applied) [#389](https://github.com/nf-core/rnafusion/pull/389)
- Now publishing convert2bed output to convert2bed to keep the output file [#420](https://github.com/nf-core/rnafusion/pull/420)
- No more checks for existence of samplesheet, which made building references fail (building references uses a fake sample sheet if none is provided) [#420](https://github.com/nf-core/rnafusion/pull/420)
- `--annotate --examine_coding_effect` to collect more data from fusioninspector [#426](https://github.com/nf-core/rnafusion/pull/426)
- Update vcf creation to get positions/chromosomes and strands even when fusions are filtered out by fusioninspector, using the csv output from fusion-report [#443](https://github.com/nf-core/rnafusion/pull/443)
- `Arriba` updated to 2.4.0 [#429](https://github.com/nf-core/rnafusion/pull/429)
- Change megafusion into vcf_collect, taking into account e.g. the annotation and coding effects outputs from fusioninspector, HGNC ids, frame status... [#414](https://github.com/nf-core/rnafusion/pull/414)
- CI tests on `--all` instead of each tool separately, and include trimmed/not trimmed matrix tests [#430](https://github.com/nf-core/rnafusion/pull/430)
- AWS tests on `--all` instead of each tool separately, and include trimmed/not trimmed matrix tests [#433](https://github.com/nf-core/rnafusion/pull/433)
- Update `fusion-report` to 2.1.8, updated COSMIC database to fix 404 error, fix download of references via proxy and removing FusionGDB database [#445](https://github.com/nf-core/rnafusion/pull/445)
- Update documentation [#446](https://github.com/nf-core/rnafusion/pull/446)

### Fixed

- Fix channel i/o issue in StringTie workflow and add StringTie in github CI tests [#416](https://github.com/nf-core/rnafusion/pull/416)
- Update modules, and make sure MultiQC displays the QC results properly [#440](https://github.com/nf-core/rnafusion/pull/440)
- Add 'when' condition to run CollectInsertSizeMetrics only when STAR-fusion bam files are available [#444](https://github.com/nf-core/rnafusion/pull/444)

### Removed

- Remove `squid` and `pizzly` fusion detection tools [#406](https://github.com/nf-core/rnafusion/pull/406)
- Remove harsh trimming option `--trim` [#413](https://github.com/nf-core/rnafusion/pull/413)
- Remove qualimap rna_seq [#407](https://github.com/nf-core/rnafusion/pull/407)

## v2.4.0 - [2023/09/22]

### Added

### Changed

- Use institutional configs by default [#381](https://github.com/nf-core/rnafusion/pull/381)
- Remove redundant indexing in starfusion and qc workflows [#387](https://github.com/nf-core/rnafusion/pull/387)
- Output bai files in same directory as bam files [#387](https://github.com/nf-core/rnafusion/pull/387)
- Update and review documentation [#396](https://github.com/nf-core/rnafusion/pull/396)
- Update picard container for `PICARD_COLLECTRNASEQMETRICS` to 3.0.0 [#395](https://github.com/nf-core/rnafusion/pull/395)
- Renamed output files [#395](https://github.com/nf-core/rnafusion/pull/395)
  - `Arriba` visualisation pdf from meta.id to meta.id_combined_fusions_arriba_visualisation
  - cram file from output bam of `STAR_FOR_ARRIBA`: meta.id to meta.id_star_for_arriba
  - cram file from output bam of `STAR_FOR_STARFUSION`: meta.id to meta.id.star_for_starfusion.Aligned.sortedByCoord.out
  - `fusion-report` index.html file to meta.id_fusionreport_index.html
  - meta.id.vcf output from `MEGAFUSION` to meta.id_fusion_data.vcf
  - Update metro map [#428](https://github.com/nf-core/rnafusion/pull/428)

### Fixed

- Tail trimming for reverse reads [#379](https://github.com/nf-core/rnafusion/pull/379)
- Set html files as optional in fusionreport [#380](https://github.com/nf-core/rnafusion/pull/380)
- Provide gene count file by default when running STAR_FOR_STARFUSION [#385](https://github.com/nf-core/rnafusion/pull/385)
- Fix fusion-report issue with MACOXS directories [#386](https://github.com/nf-core/rnafusion/pull/386)
- The fusion lists is updated to contain two branches, one in case no fusions are detected and one for if fusions are detected, that will be used to feed to fusioninspector, megafusion, arriba visualisation [#388](https://github.com/nf-core/rnafusion/pull/388)
- Update fusionreport to 2.1.5p4 to fix 403 error in downloading databases [#403](https://github.com/nf-core/rnafusion/pull/403)

### Removed

- `samtools sort` and `samtools index` for `arriba` workflow were dispensable and were removed [#395](https://github.com/nf-core/rnafusion/pull/395)
- Removed trimmed fastqc report from multiqc [#394](https://github.com/nf-core/rnafusion/pull/394)

## v2.3.0 - [2023/04/24]

### Added

- Shell specification to bash
- COSMIC password put into quotes
- Trimmed reads QC in MultiQC
- Add `ARRIBA_VISUALISATION` to processed affected by `--skip_vis`
- Option `fusionreport_filter` to in/activate fusionreport displaying of fusions detected by 2 or more tools

### Changed

- `Arriba` visualisation now runs for FusionInspector (combined tools) results, not only `Arriba` results
- Updated metro map with trimming options and placed `Arriba` visualisation after `FusionInspector`
- Exit with error when using squid in combination with any ensembl version different from 102

### Fixed

- Channel issue with indexing of files with using `--cram squid`
- `Arriba` references published in the correct folder

### Removed

## v2.2.0 - [2023/03/13]

### Added

- exitStatus 140 now part of the retry strategy
- stubs to all local modules
- `--stringtie` option added with StringTie v2.2.1 to detect splicing events. Not included in `fusion-report` or `fusionInspector` summaries. Included in the `--all` workflow
- Generation of ribosomal RNA interval list with build_references and use it in picard CollectRnaMetrics
- Add csv output to fusionreport
- Trimming workflow using `fastp`: use trimmed reads for all tools
- `whitelist` parameter to add custom fusions to the detected ones and consider the whole for the `fusionInspector` analysis
- Compression to CRAM files for arriba, squid and starfusion workflows (fusioncatcher and pizzly do not produce SAM/BAM files, fusioninspector BAM files are too small to benefit from compression)
- `--qiagen` option to download from QIAGEN instead of COSMIC (use QIAGEN user and password for `cosmic_username` and `cosmic_passwd`)
- Bumped `STAR genomegenerate` time request for building as it was always crashing for most users
- Fixed issue with arriba visualisation parameters [#326](https://github.com/nf-core/rnafusion/issues/326)

### Changed

- Test profiles unified under 'test' but if the references do not all need to be downloaded, run with `-stub`
- Update CUSTOM_DUMPSOFTWAREVERSIONS to use multiqc version 1.13
- Updated to nf-core template 2.7.2, with all module updates
- `MultiQC` updated to 1.13a in process dumpsoftwareversion
- Patch fusion-report version with fixed mittelman DB and DB extraction date written into software_version.yaml
- `Arriba` references back to downloading with `build_references` instead of taking from container
- `Arriba` visualisation now running with `Arriba` v2.3.0
- Updated `STAR-Fusion` to 1.12.0

### Fixed

- AWS megatest to display on nf-core website
- `arriba` visualisation references updated to 2.3.0
- Removed issue with multiple outputs in samtools view for squid

### Removed

- FUSIONINSPECTOR_DEV process as the option fusioninspector_limitSjdbInsertNsj is part of the main starfusion release

## v2.1.0 - [2022/07/12]

### Added

- `FusionCatcher` single_end support for single reads ABOVE 130 bp
- `--fusioninspector_only` parameter to run FusionInspector standalone feeding gene list manually with parameter `--fusioninspector_fusions PATH`
- `--fusioncatcher_limitSjdbInsertNsj` parameter to feed --limitSjdbInsertNsj to FusionCatcher
- `--fusioninspector_limitSjdbInsertNsj` parameter to feed --limitSjdbInsertNsj to FusionInspector !!Any other value than default will use the dev version of FusionInspector!!
- OPTIONAL trimming option `--trim` for hard trimming to 75 bp in case of high read-through. Only fusioncatcher uses trimmed reads as STAR-based fusion detection tools are less sensitive to read-through
- `picard` metrics, STAR final log, and QualiMap output included in `MultiQC` report

### Changed

- `seq_platform` and `seq_center` changed from boolean to string
- `seq_platform` set to an empty string and `seq_center` set to an empty string if not existing
- Arriba use ensembl references-built starindex independently of `starfusion_build` parameter
- ftp to http protocol for STARFUSION_BUILD process `Pfam-A.hmm.gz` download as ftp causes issues on some servers
- Updated README and usage documentation with more detailed information and metro map
- Arriba use ensembl references-built starindex independently of starfusion_build parameter
- Update of the single-end reads support table in README, added recommendation to use single-end reads only in last resort
- STAR updated to 2.7.10a
- Arriba updated to 2.3.0, references for blacklist and protein domains changed to 2.3.0 from singularity/docker container -> arriba download of references not necessary any more
- multiQC updated to 1.13a
- picard updated to 2.27.4
- dumpsoftwareversions module updated to use multiqc=1.12 containers

### Fixed

- FusionInspector does not mix sample reads with fusion lists and meta information from other samples anymore
- Arriba visualisation does not mix sample reads with fusion lists and meta information from other samples anymore
- logging of STAR-fusion and fusionreport version

### Removed

## v2.0.0 - [2022/05/19]

Update to DSL2 and newer software/reference versions

### Added

- Added `qualimap/rnaseq v2.2.2d` from nf-core modules
- Added UCSC `gtfToGenePred v377`
- Added `picard CollectRnaSeqMetrics v2.26.10`
- Added `picard MarkDuplicates v2.26.10` from nf-core modules
- Added `cat/fastqc` from nf-core modules
- Added possibility for manually feeding the results of fusions from different tools to speed-up reruns
- STAR-Fusion references can be downloaded or built but downloaded references are NOT RECOMMENDED as not thoroughly tested (--starfusion_build parameter is true by default, use --starfusion_build false to use downloaded STAR-Fusion references).

### Changed

- Upgrade default ensembl version to `102`
- Upgrade to `nf-core/tools v2.3.2`
- Upgrade `Arriba v1.2.0` to `Arriba v2.2.1`
- Upgrade `FusionCatcher v1.20` to `FusionCatcher v1.33`
- Upgrade `STAR-fusion v1.8.1` to `STAR-fusion v1.10.1`
- Upgrade `STAR v2.7.1` to `STAR v2.7.9`
- Upgrade `fusion-report v2.1.3` to `fusion-report v2.1.5`
- Upgrade `kallisto v0.44.0` to `kallisto v0.46.2`
- Upgrade `fastqc v0.11.8` to `fastqc v0.11.9`
- Upgrade `samtools v1.9` to `samtools v1.15.1`
- Upgrade `arriba` references from `v1.2.0` to `v2.1.0`
- Upgrade `fusioncatcher` references from `v98` to `v102`
- Use `arriba` (detect only), `kallisto` and `STAR` from nf-core modules
- Instead of separate script to build the references, added `--build_references` argument in the main
- `--fasta` argument is not required with `--build_references` and set by default to the ensembl references built in the detection workflow
- CI test done on stubs of reference building for subprocesses ensembl and arriba

Parameters for `STAR` for `arriba` changed from:

```bash
--readFilesCommand zcat \\
        --outSAMtype BAM Unsorted \\
--outStd BAM_Unsorted \\
--outSAMunmapped Within \\
--outBAMcompression 0 \\
--outFilterMultimapNmax 1 \\
--outFilterMismatchNmax 3 \\
--chimSegmentMin 10 \\
--chimOutType WithinBAM SoftClip \\
--chimJunctionOverhangMin 10 \\
--chimScoreMin 1 \\
--chimScoreDropMax 30 \\
--chimScoreJunctionNonGTAG 0 \\
--chimScoreSeparation 1 \\
--alignSJstitchMismatchNmax 5 -1 5 5 \\
--chimSegmentReadGapMax 3 \\
--sjdbOverhang ${params.read_length - 1}
```

to

```bash
--readFilesCommand zcat \
--outSAMtype BAM Unsorted \
--outSAMunmapped Within \
--outBAMcompression 0 \
--outFilterMultimapNmax 50 \
--peOverlapNbasesMin 10 \
--alignSplicedMateMapLminOverLmate 0.5 \
--alignSJstitchMismatchNmax 5 -1 5 5 \
--chimSegmentMin 10 \
--chimOutType WithinBAM HardClip \
--chimJunctionOverhangMin 10 \
--chimScoreDropMax 30 \
--chimScoreJunctionNonGTAG 0 \
--chimScoreSeparation 1 \
--chimSegmentReadGapMax 3 \
--chimMultimapNmax 50
```

As recommended [here](https://arriba.readthedocs.io/en/latest/workflow/).

Parameters for `STAR` for `STAR-fusion` changed from:

```bash
--twopassMode Basic \\
--outReadsUnmapped None \\
--chimSegmentMin 12 \\
--chimJunctionOverhangMin 12 \\
--alignSJDBoverhangMin 10 \\
--alignMatesGapMax 100000 \\
--alignIntronMax 100000 \\
--chimSegmentReadGapMax 3 \\
--alignSJstitchMismatchNmax 5 -1 5 5 \\
--runThreadN ${task.cpus} \\
--outSAMstrandField intronMotif ${avail_mem} \\
--outSAMunmapped Within \\
--outSAMtype BAM Unsorted \\
--outSAMattrRGline ID:GRPundef \\
--chimMultimapScoreRange 10 \\
--chimMultimapNmax 10 \\
--chimNonchimScoreDropMin 10 \\
--peOverlapNbasesMin 12 \\
--peOverlapMMp 0.1 \\
--readFilesCommand zcat \\
--sjdbOverhang ${params.read_length - 1} \\
--chimOutJunctionFormat 1
```

to

```bash
--outReadsUnmapped None \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outSAMstrandField intronMotif \
--outSAMunmapped Within \
--chimSegmentMin 12 \
--chimJunctionOverhangMin 8 \
--chimOutJunctionFormat 1 \
--alignSJDBoverhangMin 10 \
--alignMatesGapMax 100000 \
--alignIntronMax 100000 \
--alignSJstitchMismatchNmax 5 -1 5 5 \
--chimMultimapScoreRange 3 \
--chimScoreJunctionNonGTAG -4 \
--chimMultimapNmax 20 \
--chimNonchimScoreDropMin 10 \
--peOverlapNbasesMin 12 \
--peOverlapMMp 0.1 \
--alignInsertionFlush Right \
--alignSplicedMateMapLminOverLmate 0 \
--alignSplicedMateMapLmin 30 \
--chimOutType Junctions
```

`Homo_sapiens.${params.genome}.${ensembl_version}.gtf.gz` used for squid and arriba, `Homo_sapiens.${params.genome}.${ensembl_version}.chr.gtf.gz` used for STAR-fusion and the quality control as the quality control is based on the STAR-fusion alignment.

### Fixed

### Removed

- Ericscript tool
- GRCh37 support. Subdirectory with params.genome are removed
- Running with conda

## v1.3.0 - [2020/07/15]

- Using official STAR-Fusion container [#160](https://github.com/nf-core/rnafusion/issues/160)

### Added

- Added social preview image [#107](https://github.com/nf-core/rnafusion/issues/107)
- Added support for GRCh37 genome assembly [#77](https://github.com/nf-core/rnafusion/issues/77)

### Changed

- Upgrade `fusion-report v2.1.2` to `fusion-report v2.1.3`
- Upgrade `fusion-report v2.1.1` to `fusion-report v2.1.2`
- Upgrade `fusion-report v2.1.0` to `fusion-report v2.1.1`
- Upgrade `Arriba v1.1.0` to `Arriba v1.2.0`
- Upgrade `fusion-report v2.0.2` to `fusion-report v2.1.0`

### Fixed

- Missing `strip-components` in `download-references.nf/star-fusion` [#148](https://github.com/nf-core/rnafusion/issues/148)
- Missing version prefix for cdna [#143](https://github.com/nf-core/rnafusion/issues/143)
- `samtools` missing header in empty file for FusionInspector [ref](https://github.com/STAR-Fusion/STAR-Fusion/issues/191)
- Removed `profile` from helper scripts [#139](https://github.com/nf-core/rnafusion/issues/139)
- Wrong url path for `Pfam-A.hmm.gz` [#140](https://github.com/nf-core/rnafusion/issues/140)

### Removed

- Removed `scripts/download-singularity-img.sh` and `download-singularity-img.nf` as they are not necessary any more

---

## v1.1.0 - [2020/02/10]

- Fusion gene detection tools:
  - `Arriba v1.1.0`
  - `Ericscript v0.5.5`
  - `Fusioncatcher v1.20`
  - `Pizzly v0.37.3`
  - `Squid v1.5`
  - `STAR-Fusion v1.6.0`
- Visualization tools:
  - `Arriba v1.1.0`
  - `FusionInspector v1.3.1`
- Other tools:
  - `fusion-report v2.0.1`
  - `FastQ v0.11.8`
  - `MultiQC v1.7`
  - `STAR aligner v2.7.0f`

### Added

- Added `Arriba 1.1.0` [#63](https://github.com/nf-core/rnafusion/issues/63)
- Added Batch mode [#54](https://github.com/nf-core/rnafusion/issues/54)

### Changed

- Updated examples and configurations
- Upgraded `fusion-report v1.0.0` to `fusion-report v2.0.1`
- Divided `running_tools` into fusion and visualization tools
- Updated `STAR` in `Squid`, `Fusion-Inspector` version to `2.7.0f`
- Upgraded `STAR-Fusion v1.5.0` to `STAR-Fusion v1.6.0` [#83](https://github.com/nf-core/rnafusion/issues/83)
- Parameter `igenomesIgnore` renamed to `igenome` [#81](https://github.com/nf-core/rnafusion/issues/81)
- Finished STAR-Fusion file renaming [#18](https://github.com/nf-core/rnafusion/issues/18)
- Updated logos
- Updated to nf-core `1.8` TEMPLATE

### Fixed

- iGenomes optional, but not really [#91](https://github.com/nf-core/rnafusion/issues/91)
- Updated `fusioncatcher` to latest `1.20` version also solving [#95](https://github.com/nf-core/rnafusion/issues/95)

### Removed

- Variables `pizzly_fasta` and `pizzly_gtf` have been removed and replaced with `transcript` and `gtf`
- `Jenkisfile`, test configuration, pylintrc configuration
- Removed `igenomes.config` because the pipeline only supports `Ensembl` version

---

## v1.0.2 - [2019/05/13]

### Changed

- Bumped nf-core template to 1.6 [#69](https://github.com/nf-core/rnafusion/pull/69)

### Fixed

- Fixed COSMIC parameters not wrapped in quotes [#75](https://github.com/nf-core/rnafusion/issues/75)
- Implemented output output for fusion tools [#72](https://github.com/nf-core/rnafusion/issues/72)
- Fixed reference download link for STAR-Fusion [#71](https://github.com/nf-core/rnafusion/issues/71)

---

## v1.0.1 - [2019/04/06]

### Added

- Added support for extra parameters for tools STAR-Fusion, FusionCatcher and fusion-report
- Added example configuration for `singularity` and `docker`
- Added [fusion-report](https://github.com/matq007/fusion-report) into the stack [#62](https://github.com/nf-core/rnafusion/issues/62), [#55](https://github.com/nf-core/rnafusion/issues/55), [#53](https://github.com/nf-core/rnafusion/issues/53), [#51](https://github.com/nf-core/rnafusion/issues/51)
- Added nextflow helper script `download-singularity-img.nf`
- Added nextflow helper script `download-references.nf`
- Added `Jenkinsfile` for in-house testing

### Changed

- Updated installation of `FusionCatcher` (available now on bioconda)

### Fixed

- Fixed empty symlinks (`input.X`) in fusion-report [#68](https://github.com/nf-core/rnafusion/issues/68)
- Fixed FASTA issues [#60](https://github.com/nf-core/rnafusion/issues/60)
- Fixed centralized nf-core/config [#64](https://github.com/nf-core/rnafusion/issues/64)
- Fixed `scrape_software_versions.py` to parse tools versions correctly [#65](https://github.com/nf-core/rnafusion/issues/65)

### Removed

- Removed `Singularity`

---

## v1.0 - [2018/02/14]

Version 1.0 marks the first production release of this pipeline under the nf-core flag.
The pipeline includes additional help scripts to download references for fusion tools and Singularity images.

Initial release of nf-core/rnafusion, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
