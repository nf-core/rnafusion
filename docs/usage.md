# nf-core/rnafusion: Usage <!-- omit in toc -->

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/rnafusion/usage](https://nf-co.re/rnafusion/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

- [:warning: Please read this documentation on the nf-core website: https://nf-co.re/rnafusion/usage](#warning-please-read-this-documentation-on-the-nf-core-website-httpsnf-corernafusionusage)
- [Introduction](#introduction)
- [Download references](#download-references)
    - [Download all references](#download-all-references)
    - [Download specific references](#download-specific-references)
    - [Download and build CTAT](#download-and-build-ctat)
    - [Download GRCh37 references](#download-grch37-references)
    - [Tool reference requirements](#tool-reference-requirements)
- [Running the pipeline](#running-the-pipeline)
    - [Updating the pipeline](#updating-the-pipeline)
    - [Reproducibility](#reproducibility)
- [Core Nextflow arguments](#core-nextflow-arguments)
    - [`-profile`](#-profile)
    - [`-resume`](#-resume)
    - [`-c`](#-c)
        - [Custom resource requests](#custom-resource-requests)
    - [Running in the background](#running-in-the-background)
        - [Nextflow memory requirements](#nextflow-memory-requirements)
- [Pipeline specific arguments](#pipeline-specific-arguments)
    - [--input](#--input)
    - [--single_end](#--single_end)
    - [Tool flags](#tool-flags)
    - [--arriba](#--arriba)
    - [--ericscript](#--ericscript)
    - [--fusioncatcher](#--fusioncatcher)
    - [--fusion_report](#--fusion_report)
    - [--pizzly](#--pizzly)
    - [--squid](#--squid)
    - [--star_fusion](#--star_fusion)
- [Visualization flags](#visualization-flags)
    - [--arriba_vis](#--arriba_vis)
    - [--fusion_inspector](#--fusion_inspector)
- [Reference genomes](#reference-genomes)
    - [--arriba_ref](#--arriba_ref)
    - [--databases](#--databases)
    - [--ericscript_ref](#--ericscript_ref)
    - [--fasta](#--fasta)
    - [--fusioncatcher_ref](#--fusioncatcher_ref)
    - [--genome](#--genome)
    - [--gtf](#--gtf)
    - [--reference_release](#--reference_release)
    - [--star_index](#--star_index)
    - [--star_fusion_ref](#--star_fusion_ref)
    - [--transcript](#--transcript)
- [Other command line parameters](#other-command-line-parameters)
    - [--debug](#--debug)
    - [--read_length](#--read_length)
    - [--outdir](#--outdir)
    - [--email](#--email)
    - [--email_on_fail](#--email_on_fail)
    - [--max_multiqc_email_size](#--max_multiqc_email_size)
    - [-name](#-name)
    - [--custom_config_version](#--custom_config_version)
    - [--custom_config_base](#--custom_config_base)
- [Job resources](#job-resources)
    - [Automatic resubmission](#automatic-resubmission)

## Download references

The rnafusion pipeline needs references for the fusion detection tools, so downloading these is a prerequisite.

Downloading references manually is a tedious long process.
To make the pipeline easier to work with, we provide a script to download all necessary references for fusion detection tools.

> **TL;DR:** Make sure to download the correct references for your need!

### Download all references

```bash
# Replace <COSMIC_USER> and <COSMIC_PASSWD> with yout credentials from COSMIC
nextflow run nf-core/rnafusion/download-references.nf \
  --download_all \
  --outdir <PATH> \
  --cosmic_usr <COSMIC_USER> --cosmic_passwd <COSMIC_PASSWD>
```

### Download specific references

```bash
# Example of downloading base references

nextflow run nf-core/rnafusion/download-references.nf \
--base \
--outdir <PATH>

# Example of downloading references for specific tool

nextflow run nf-core/rnafusion/download-references.nf \
--arriba \
--outdir <PATH>

nextflow run nf-core/rnafusion/download-references.nf \
--star_fusion \
--outdir <PATH>

nextflow run nf-core/rnafusion/download-references.nf \
--fusioncatcher \
--outdir <PATH>

nextflow run nf-core/rnafusion/download-references.nf \
--ericscript \
--outdir <PATH>

nextflow run nf-core/rnafusion/download-references.nf \
--fusion_report \
--outdir <PATH>
  --cosmic_usr <COSMIC_USER> --cosmic_passwd <COSMIC_PASSWD>
```

### Download and build CTAT

```bash
nextflow run nf-core/rnafusion/build-ctat.nf \
  --genome GRCh38 \
  --outdir <PATH> \
  --fasta <PATH>/<FASTA \
  --gtf <PATH>/<GTF>
```

### Download GRCh37 references

```bash
# GRCh38 genome assembly is used by default.
# To use the previous assembly specify it using the --genome flag
nextflow run nf-core/rnafusion/download-references.nf \
  --genome GRCh37 \
  --download_all \
  --outdir <PATH> \
  --cosmic_usr <COSMIC_USER> --cosmic_passwd <COSMIC_PASSWD>

# Please note that using the above example command downloads NCBI-based references for STAR-Fusion.
# To use Ensembl-based references run the following command with the same <PATH> as used above

nextflow run nf-core/rnafusion/build-ctat.nf \
  --genome GRCh37 \
  --outdir <PATH> \
  --fasta <PATH>/<FASTA> \
  --gtf <PATH>/<GTF>
```

### Tool reference requirements

| Tool             |        FASTA       |         GTF        |     STAR-index     |     Genome     |       Other       |
| ---------------- | :----------------: | :----------------: | :----------------: | :------------: | :----------------: |
| Arriba           | :white_check_mark: | :white_check_mark: | :white_check_mark: | GRCh37, GRCh38 | `custom_reference` |
| EricScript       |         :x:        |         :x:        |         :x:        | GRCh37, GRCh38 | `custom_reference` |
| FusionCatcher    |         :x:        |         :x:        |         :x:        |     GRCh38     | `custom_reference` |
| Fusion-Inspector | :white_check_mark: | :white_check_mark: | :white_check_mark: | GRCh37, GRCh38 |  `ctat_genome_lib` |
| fusion-report    |         :x:        |         :x:        |         :x:        | GRCh37, GRCh38 |     `databases`    |
| Pizzly           |         :x:        | :white_check_mark: | :white_check_mark: | GRCh37, GRCh38 |       `cDNA`       |
| Squid            |         :x:        | :white_check_mark: | :white_check_mark: | GRCh37, GRCh38 |          -         |
| Star-Fusion      | :white_check_mark: | :white_check_mark: | :white_check_mark: | GRCh37, GRCh38 |  `ctat_genome_lib` |

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```console
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```console
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 3 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice.

```console
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz,
```

| Column         | Description                                                                                                                                                                            |
|----------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `sample`       | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1`      | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2`      | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows.

```console
nextflow run nf-core/rnafusion --input samplesheet.csv --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```console
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```console
nextflow pull nf-core/rnafusion
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/rnafusion releases page](https://github.com/nf-core/rnafusion/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below. When using Biocontainers, most of these software packaging methods pull Docker containers from quay.io e.g [FastQC](https://quay.io/repository/biocontainers/fastqc) except for Singularity which directly downloads Singularity images via https hosted by the [Galaxy project](https://depot.galaxyproject.org/singularity/) and Conda which downloads and installs software locally from [Bioconda](https://bioconda.github.io/).

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

- `docker`
    - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
    - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
    - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
    - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
    - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `conda`
    - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
- `test`
    - A profile with a complete configuration for automated testing
    - Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/rnaseq pipeline is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

Command executed:
    STAR \
        --genomeDir star \
        --readFilesIn WT_REP1_trimmed.fq.gz  \
        --runThreadN 2 \
        --outFileNamePrefix WT_REP1. \
        <TRUNCATED>

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    STAR --genomeDir star --readFilesIn WT_REP1_trimmed.fq.gz --runThreadN 2 --outFileNamePrefix WT_REP1. <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process. The quickest way is to search for `process STAR_ALIGN` in the [nf-core/rnaseq Github repo](https://github.com/nf-core/rnaseq/search?q=process+STAR_ALIGN). We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so based on the search results the file we want is `modules/nf-core/software/star/align/main.nf`. If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_high`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L9). The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements. The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L33-L37) which in this case is defined as 72GB. Providing you haven't set any other standard nf-core parameters to __cap__ the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `STAR_ALIGN` process failure by creating a custom config file that sets at least 72GB of memory, in this case increased to 100GB. The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: STAR_ALIGN {
        memory = 100.GB
    }
}
```

> **NB:** We specify just the process name i.e. `STAR_ALIGN` in the config file and not the full task name string that is printed to screen in the error message or on the terminal whilst the pipeline is running i.e. `RNASEQ:ALIGN_STAR:STAR_ALIGN`. You may get a warning suggesting that the process selector isn't recognised but you can ignore that if the process name has been specified correctly. This is something that needs to be fixed upstream in core Nextflow.

### Updating containers

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

    - For Docker:

        ```nextflow
        process {
            withName: PANGOLIN {
                container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
            }
        }
        ```

    - For Singularity:

        ```nextflow
        process {
            withName: PANGOLIN {
                container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
            }
        }
        ```

    - For Conda:

        ```nextflow
        process {
            withName: PANGOLIN {
                conda = 'bioconda::pangolin=3.0.5'
            }
        }
        ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```console
NXF_OPTS='-Xms1g -Xmx4g'
```

## Pipeline specific arguments

### --input

Use this to specify the location of your input FastQ files. For example:

```bash
--input 'path/to/data/sample_*_{1,2}.fastq.gz'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

### --single_end

By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--single_end` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--input`. For example:

```bash
--single_end --input '*.fastq'
```

### Tool flags

### --arriba

If enabled, executes `Arriba` tool.

- `--arriba_opt`
    - Specify additional parameters. For more information, please refer to the [documentation](http://arriba.readthedocs.io/en/latest/quickstart/) of the tool.

### --ericscript

If enabled, executes `Ericscript` tool.

- `--ericscript_opt`
    - Specify additional parameters. For more information, please refer to the [documentation](https://sites.google.com/site/bioericscript/home) of the tool.

### --fusioncatcher

If enabled, executes `Fusioncatcher` tool.

> N.B. that Fusioncatcher is not available when using the `GRCh37` genome assembly.

- `--fusioncatcher_opt`
    - Specify additional parameters. For more information, please refer to the [documentation](https://github.com/ndaniel/fusioncatcher/blob/master/doc/manual.md) of the tool.

### --fusion_report

If enabled, download databases for `fusion-report`.

- `fusion_report_opt`
    - Specify additional parameters. For more information, please refer to the [documentation](https://matq007.github.io/fusion-report/#/) of the tool.

### --pizzly

If enabled, executes `Pizzly` tool.

- `--pizzly_k`
    - Number of k-mers. Default `31`.

### --squid

If enabled, executes `Squid` tool.

### --star_fusion

If enabled, executes `STAR-Fusion` tool.

- `--star_fusion_opt`
    - Parameter for specifying additional parameters. For more information, please refer to the [documentation](https://github.com/STAR-Fusion/STAR-Fusion/wiki) of the tool.

## Visualization flags

### --arriba_vis

If enabled, executes build in `Arriba` visualization tool.

### --fusion_inspector

If enabled, executes `Fusion-Inspector` tool.

## Reference genomes

### --arriba_ref

```bash
--arriba_ref '<path to Arriba reference>'
```

### --databases

Required databases in order to run `fusion-report`.

```bash
--databases '<path to fusion-report databases>'
```

### --ericscript_ref

Required reference in order to run `EricScript`.

```bash
--ericscript_ref '<path to EricScript reference>'
```

### --fasta

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--fasta '<path to Fasta reference>'
```

### --fusioncatcher_ref

Required reference in order to run `Fusioncatcher`.

```bash
--fusioncatcher_ref '<path to Fusioncatcher reference>'
```

### --genome

This pipeline uses `Homo Sapiens` version `GRCh38` by default. Assembly `GRCh37` is optionally available.

> N.B. that using `GRCh37` precludes use of the `Fusioncatcher` tool.
> Also make sure to specify `--genomes_base`.

```bash
--genome 'GRCh38' --genome_base '</path/to/references>'
```

### --gtf

Required annotation file.

```bash
--gtf '<path to GTF annotation>'
```

### --reference_release

Ensembl version.

```bash
# ftp://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/
--reference_release '97'
```

### --star_index

If you prefer, you can specify the full path for `STAR` index when you run the pipeline.
If not specified, the pipeline will build the index using for reads with length `100bp` (can be adjusted with parameter `--read_length`).

```bash
--star_index '<path to STAR index>'
```

### --star_fusion_ref

Required reference in order to run `STAR-Fusion`.

```bash
--star_fusion_ref '<path to STAR-Fusion reference>'
```

### --transcript

Required transcript file.

```bash
--transcript '<path to transcript reference>'
```

## Other command line parameters

### --debug

To run only a specific tool (testing freshly implemented tool) just add `--debug` parameter.
This parameter only works on **fusion tools only**!

### --read_length

Length is used to build a STAR index.
Default is `100bp` (Illumina).

### --outdir

The output directory where the results will be saved.

### --email

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits.
If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### --email_on_fail

This works exactly as with `--email`, except emails are only sent if the workflow is not successful.

### --max_multiqc_email_size

Threshold size for MultiQC report to be attached in notification email.
If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB).

### -name

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### --custom_config_version

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default: `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### --custom_config_base

If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet. If you don't need them, then this is not a problem. If you do need them,
you should download the files from the repo and tell nextflow where to find them with the
`custom_config_base` option. For example:

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time.
For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original).
If it still fails after three times then the pipeline is stopped.
