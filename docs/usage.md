# nf-core/rnafusion: Usage

## Table of contents

<!-- Install Atom plugin markdown-toc-auto for this ToC to auto-update on save -->
<!-- TOC START min:2 max:3 link:true asterisk:true update:true -->
* [Table of contents](#table-of-contents)
* [Introduction](#introduction)
* [Running the pipeline](#running-the-pipeline)
  * [Using Docker](#running-the-pipeline-using-docker)
  * [Using Singularity](#running-the-pipeline-using-singularity)
* [Updating the pipeline](#updating-the-pipeline)
* [Reproducibility](#reproducibility)
* [Main arguments](#main-arguments)
  * [`-profile`](#-profile-single-dash)
    * [`awsbatch`](#awsbatch)
    * [`conda`](#conda)
    * [`docker`](#docker)
    * [`singularity`](#singularity)
    * [`test`](#test)
  * [`--reads`](#--reads)
  * [`--singleEnd`](#--singleend)
* [Tool flags](#tool-flags)
  * [`--star_fusion`](#--star_fusion)
    * [`--star_fusion_opt`](#--star_fusion_opt)
  * [`--fusioncatcher`](#--fusioncatcher)
    * [`--fusioncatcher_opt`](#--fusioncatcher_opt)
  * [`--ericscript`](#--ericscript)
  * [`--pizzly`](#--pizzly)  
  * [`--squid`](#--squid)
  * [`--fusion_report_opt`](#--fusion_report_opt)
  * [`--debug`](#--debug)
* [Visualization flags](#visualization-flags)
  * [`--fusion_inspector`](#--fusion_inspector)
* [Reference genomes](#reference-genomes)
  * [`--fasta`](#--fasta)
  * [`--gtf`](#--gtf)
  * [`--star_index`](#--star_index)
  * [`--star_fusion_ref`](#--star_fusion_ref)
  * [`--fusioncatcher_ref`](#--fusioncatcher_ref)
  * [`--ericscript_ref`](#--ericscript_ref)
  * [`--pizzly_fasta`](#--pizzly_fasta)
  * [`--pizzly_gtf`](#--pizzly_gtf)
  * [`--genome` (using iGenomes)](#--genome-using-igenomes)
  * [`--igenomesIgnore`](#--igenomesignore)
* [Job resources](#job-resources)
  * [Automatic resubmission](#automatic-resubmission)
  * [Custom resource requests](#custom-resource-requests)
* [AWS Batch specific parameters](#aws-batch-specific-parameters)
  * [`--awsqueue`](#--awsqueue)
  * [`--awsregion`](#--awsregion)
* [Other command line parameters](#other-command-line-parameters)
  * [`--read_length`](#--read_length)
  * [`--outdir`](#--outdir)
  * [`--email`](#--email)
  * [`-name`](#-name)
  * [`-resume`](#-resume)
  * [`-c`](#-c)
  * [`--custom_config_version`](#--custom_config_version)
  * [`--custom_config_base`](#--custom_config_base)
  * [`--max_memory`](#--max_memory)
  * [`--max_time`](#--max_time)
  * [`--max_cpus`](#--max_cpus)
  * [`--plaintext_email`](#--plaintext_email)
  * [`--monochrome_logs`](#--monochrome_logs)
  * [`--multiqc_config`](#--multiqc_config)
<!-- TOC END -->

## Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline

The typical command for running the whole pipeline is as follows:

### Running the pipeline using Docker

This will launch the pipeline using `docker` with configuration profile [example-docker.config](https://github.com/nf-core/rnafusion/blob/master/example/custom-docker.config). See below for more information about profiles.

```bash
# With custom fasta and gtf (Ensembl example)
nextflow run nf-core/rnafusion
  --reads '*_R{1,2}.fastq.gz'
  -profile docker -c 'example/custom-docker.config'
  --fasta 'Homo_sapiens.GRCh38.95.all.fa'
  --gtf '/Homo_sapiens.GRCh38.95.chr.gtf'
  --star_fusion
  --fusioncatcher
  --ericscript
  --pizzly
  --squid
  --fusion_inspector

# With NCBI GRCh38 genome reference
nextflow run nf-core/rnafusion
  --reads '*_R{1,2}.fastq.gz'
  -profile docker -c 'example/custom-docker.config'
  --genome GRCh38
  --star_fusion
  --fusioncatcher
  --ericscript
  --pizzly
  --squid
  --fusion_inspector
```

### Running the pipeline using Singularity

First start by downloading singularity images. Sometimes the pipeline can crash if you are not using downloaded images (might be some network issues).

```bash
nextflow run nf-core/rnafusion/download-singularity-img.nf --download_all --outdir /path

# or

cd utils && sh download-singularity-img.sh /path/to/images
```

This will launch the pipeline using `singularity` with configuration profile [example-singularity.config](https://github.com/nf-core/rnafusion/blob/master/example/custom-singularity.config). See below for more information about profiles.

```bash
# With custom fasta and gtf (Ensembl example)
nextflow run nf-core/rnafusion
  --reads '*_R{1,2}.fastq.gz'
  -profile singularity -c 'example/custom-singularity.config'
  --fasta 'Homo_sapiens.GRCh38.95.all.fa'
  --gtf '/Homo_sapiens.GRCh38.95.chr.gtf'
  --star_fusion
  --fusioncatcher
  --ericscript
  --pizzly
  --squid
  --fusion_inspector

# With NCBI GRCh38 genome reference
nextflow run nf-core/rnafusion
  --reads '*_R{1,2}.fastq.gz'
  -profile singularity -c 'example/custom-singularity.config'
  --genome GRCh38
  --star_fusion
  --fusioncatcher
  --ericscript
  --pizzly
  --squid
  --fusion_inspector
```

---

It is also possible to execute **only** specific tools:

```bash
nextflow run nf-core/rnafusion
  --reads '*_R{1,2}.fastq.gz'
  --genome GRCh38 -profile docker -c 'example/custom-docker.config'
  --fusioncatcher
  --ericscript
```

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/rnafusion
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/rnafusion releases page](https://github.com/nf-core/rnafusion/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Main arguments

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile docker` - the order of arguments is important!

If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

* `awsbatch`
  * A generic configuration profile to be used with AWS Batch.
* `conda`
  * A generic configuration profile to be used with [conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from dockerhub: [`nfcore/rnafusion`](http://hub.docker.com/r/nfcore/rnafusion/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from DockerHub: [`nfcore/rnafusion`](http://hub.docker.com/r/nfcore/rnafusion/)
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

### `--reads`

Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_{1,2}.fastq.gz'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

It is not possible to run a mixture of single-end and paired-end files in one run.

### `--singleEnd`

By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--singleEnd` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--singleEnd --reads '*.fastq.gz'
```

## Tool flags

### `--star_fusion`

If enabled, executes `STAR-Fusion` tool.

* `--star_fusion_opt`
  * Parameter for specifying additional parameters. For more info, please refer to the [documentation](https://github.com/STAR-Fusion/STAR-Fusion/wiki) of the tool.
  * **Has to be specified in custom configuration file. Will not work from a command line.**

### `--fusioncatcher`

If enabled, executes `Fusioncatcher` tool.

* `--fusioncatcher_opt`
  * Parameter for specifying additional parameters. For more info, please refer to the [documentation](https://github.com/ndaniel/fusioncatcher/blob/master/doc/manual.md) of the tool.
  * **Has to be specified in custom configuration file. Will not work from a command line.**

### `--ericscript`

If enabled, executes `Ericscript` tool.

### `--pizzly`

If enabled, executes `Pizzly` tool.

### `--squid`

If enabled, executes `Squid` tool.

### `--fusion_report_opt`

* Parameter for specifying additional parameters. For more info, please refer to the fusion-report [documentation](https://matq007.github.io/fusion-report/usage.html).
* **Has to be specified in custom configuration file. Will not work from a command line.**

### `--debug`

To run only a specific tool (testing freshly implemented tool) just add `--debug` parameter. This parameter only works on **fusion tools only**!

```bash
nextflow run nf-core/rnafusion --reads '*_R{1,2}.fastq.gz' --genome GRCh38 -profile docker --star_fusion --test
```

## Visualization flags

### `--fusion_inspector`

If enabled, executes `Fusion-Inspector` tool.

## Reference genomes

The pipeline config files come bundled with paths to the illumina iGenomes reference index files. If running with docker or AWS, the configuration is set up to use the [AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/) resource.

### `--fasta`

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--fasta '[path to Fasta reference]'
```

### `--gtf`

If you prefer, you can specify the full path to your annotation when you run the pipeline:

```bash
--gtf '[path to GTF annotation]'
```

### `--star_index`

If you prefer, you can specify the full path for `STAR` index when you run the pipeline. If not specified, the pipeline will build the index using for reads with length 100bp (can be adjusted with parameter `--read_length`).

```bash
--star_index '[path to STAR index]'
```

### `--star_fusion_ref`

Required reference in order to run `STAR-Fusion`.

```bash
--star_fusion_ref '[path to STAR-Fusion reference]'
```

### `--fusioncatcher_ref`

Required reference in order to run `Fusioncatcher`.

```bash
--fusioncatcher_ref '[path to Fusioncatcher reference]'
```

### `--ericscript_ref`

Required reference in order to run `Ericscript`.

```bash
--ericscript_ref '[path to Ericscript reference]'
```

### `--pizzly_fasta`

Required reference in order to run `Pizzly`.

```bash
--pizzly_fasta '[path to Pizzly Fasta reference]'
```

### `--pizzly_gtf`

Required reference in order to run `Pizzly`.

```bash
--pizzly_gtf '[path to Pizzly GTF annotation]'
```

### `--genome` (using iGenomes)

There are 31 different species supported in the iGenomes references. To run the pipeline, you must specify which to use with the `--genome` flag.
You can find the keys to specify the genomes in the [iGenomes config file](../conf/igenomes.config). Common genomes that are supported are:

* Human
  * `--genome GRCh38` (recommended)

> **TL;DR:** The pipeline only supports Homo Sapiens. We recommend using fasta nad gtf from Ensembl database and build custom STAR-Fusion reference. Most of the tools references are based on Ensembl.

Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the iGenomes resource. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to save such a file.

The syntax for this reference configuration is as follows:

```nextflow
params {
  genomes {
    'GRCh38' {
      fasta   = '<path to the genome fasta file>' // Used if no star index given
    }
    // Any number of additional genomes, key is used with --genome
  }
}
```

### `--igenomesIgnore`

Do not load `igenomes.config` when running the pipeline. You may choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config`.

## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-core-invite.herokuapp.com/).

## AWS Batch specific parameters

Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use the `-awsbatch` profile and then specify all of the following parameters.

### `--awsqueue`

The JobQueue that you intend to use on AWS Batch.

### `--awsregion`

The AWS region to run your job in. Default is set to `eu-west-1` but can be adjusted to your needs.
Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

### `--read_length`

Length is used to build a STAR index. Default is 100bp (Illumina).

### `--outdir`

The output directory where the results will be saved.

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.
You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--custom_config_version`

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default is set to `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--custom_config_base`

If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet. If you don't need them, then this is not a problem. If you do need them,
you should download the files from the repo and tell nextflow where to find them with the
`custom_config_base` option. For example:

```bash
## Download and unzip the config files
cd /path/to/my/configs
wget https://github.com/nf-core/configs/archive/master.zip
unzip master.zip

## Run the pipeline
cd /path/to/my/data
nextflow run /path/to/pipeline/ --custom_config_base /path/to/my/configs/configs-master/
```

> Note that the nf-core/tools helper package has a `download` command to download all required pipeline
> files + singularity containers + institutional configs in one go for you, to make this process easier.

### `--max_memory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`

Set to receive plain-text e-mails instead of HTML formatted.

### `--monochrome_logs`

Set to disable colourful command line output and live life in monochrome.

### `--multiqc_config`

Specify a path to a custom MultiQC configuration file.
