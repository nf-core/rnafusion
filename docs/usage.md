<!-- omit in toc -->
# nf-core/rnafusion: Usage

<!-- omit in toc -->
## Table of contents

- [Introduction](#introduction)
- [Running the pipeline](#running-the-pipeline)
  - [Running the pipeline using Docker](#running-the-pipeline-using-docker)
  - [Running the pipeline using Singularity](#running-the-pipeline-using-singularity)
  - [Running specific tools](#running-specific-tools)
  - [Updating the pipeline](#updating-the-pipeline)
  - [Reproducibility](#reproducibility)
- [Main arguments](#main-arguments)
  - [`-profile`](#-profile)
  - [`--reads`](#--reads)
  - [`--single_end`](#--single_end)
- [Tool flags](#tool-flags)
  - [`--arriba`](#--arriba)
  - [`--ericscript`](#--ericscript)
  - [`--fusioncatcher`](#--fusioncatcher)
  - [`--fusion_report`](#--fusion_report)
  - [`--pizzly`](#--pizzly)
  - [`--squid`](#--squid)
  - [`--star_fusion`](#--star_fusion)
- [Visualization flags](#visualization-flags)
  - [`--arriba_vis`](#--arriba_vis)
  - [`--fusion_inspector`](#--fusion_inspector)
- [Reference genomes](#reference-genomes)
  - [`--arriba_ref`](#--arriba_ref)
  - [`--databases`](#--databases)
  - [`--ericscript_ref`](#--ericscript_ref)
  - [`--fasta`](#--fasta)
  - [`--fusioncatcher_ref`](#--fusioncatcher_ref)
  - [`--genome`](#--genome)
  - [`--gtf`](#--gtf)
  - [`--reference_release`](#--reference_release)
  - [`--star_index`](#--star_index)
  - [`--star_fusion_ref`](#--star_fusion_ref)
  - [`--transcript`](#--transcript)
- [Job resources](#job-resources)
  - [Automatic resubmission](#automatic-resubmission)
  - [Custom resource requests](#custom-resource-requests)
- [AWS Batch specific parameters](#aws-batch-specific-parameters)
  - [`--awsqueue`](#--awsqueue)
  - [`--awsregion`](#--awsregion)
  - [`--awscli`](#--awscli)
- [Other command line parameters](#other-command-line-parameters)
  - [`--debug`](#--debug)
  - [`--read_length`](#--read_length)
  - [`--outdir`](#--outdir)
  - [`--email`](#--email)
  - [`--email_on_fail`](#--email_on_fail)
  - [`--max_multiqc_email_size`](#--max_multiqc_email_size)
  - [`-name`](#-name)
  - [`-resume`](#-resume)
  - [`-c`](#-c)
  - [`--custom_config_version`](#--custom_config_version)
  - [`--custom_config_base`](#--custom_config_base)
  - [`--max_memory`](#--max_memory)
  - [`--max_time`](#--max_time)
  - [`--max_cpus`](#--max_cpus)
  - [`--plaintext_email`](#--plaintext_email)
  - [`--monochrome_logs`](#--monochrome_logs)
  - [`--multiqc_config`](#--multiqc_config)

## Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline

The typical command for running the pipeline is as follows.

### Running the pipeline using Docker

```bash
nextflow run nf-core/rnafusion \
  -profile docker \
  --reads '*_R{1,2}.fastq.gz' \
  --arriba \
  --star_fusion \
  --fusioncatcher \
  --ericscript \
  --pizzly \
  --squid \
  --arriba_vis \
  --fusion_inspector
```

### Running the pipeline using Singularity

```bash
nextflow run nf-core/rnafusion/download-singularity-img.nf --download_all --outdir /path
```

If the nextflow download script crashes (network issue), please use the bash script instead.

```bash
cd utils && sh download-singularity-img.sh /path/to/images
```

The command bellow will launch the pipeline using `singularity`.

```bash
nextflow run nf-core/rnafusion \
  -profile singularity \
  --reads '*_R{1,2}.fastq.gz' \
  --arriba \
  --star_fusion \
  --fusioncatcher \
  --ericscript \
  --pizzly \
  --squid \
  --arriba_vis \
  --fusion_inspector
```

### Running specific tools

```bash
nextflow run nf-core/rnafusion \
  -profile singularity -c 'example/custom-singularity.config' \
  --reads '*_R{1,2}.fastq.gz' \
  --arriba \
  --squid
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

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

- `docker`
  - A generic configuration profile to be used with [Docker](http://docker.com/)
  - Pulls software from DockerHub: [`nfcore/rnafusion`](http://hub.docker.com/r/nfcore/rnafusion/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  - Pulls software from DockerHub: [`nfcore/rnafusion`](http://hub.docker.com/r/nfcore/rnafusion/)
- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters

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

### `--single_end`

By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--single_end` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--single_end --reads '*.fastq'
```

## Tool flags

### `--arriba`

If enabled, executes `Arriba` tool.

- `--arriba_opt`
  - Specify additional parameters. For more info, please refer to the [documentation](http://arriba.readthedocs.io/en/latest/quickstart/) of the tool.

### `--ericscript`

If enabled, executes `Ericscript` tool.

- `--ericscript_opt`
  - Specify additional parameters. For more info, please refer to the [documentation](https://sites.google.com/site/bioericscript/home) of the tool.

### `--fusioncatcher`

If enabled, executes `Fusioncatcher` tool.

- `--fusioncatcher_opt`
  - Specify additional parameters. For more info, please refer to the [documentation](https://github.com/ndaniel/fusioncatcher/blob/master/doc/manual.md) of the tool.

### `--fusion_report`

If enabled, download databases for `fusion-report`.

- `fusion_report_opt`
  - Specify additional parameters. For more info, please refer to the [documentation](https://matq007.github.io/fusion-report/#/) of the tool.

### `--pizzly`

If enabled, executes `Pizzly` tool.

- `--pizzly_k`
  - Number of k-mers. Deafult 31.

### `--squid`

If enabled, executes `Squid` tool.

### `--star_fusion`

If enabled, executes `STAR-Fusion` tool.

- `--star_fusion_opt`
  - Parameter for specifying additional parameters. For more info, please refer to the [documentation](https://github.com/STAR-Fusion/STAR-Fusion/wiki) of the tool.

## Visualization flags

### `--arriba_vis`

If enabled, executes build in `Arriba` visualization tool.

### `--fusion_inspector`

If enabled, executes `Fusion-Inspector` tool.

## Reference genomes

### `--arriba_ref`

```bash
--arriba_ref '[path to Arriba reference]'
```

### `--databases`

Required databases in order to run `fusion-report`.

```bash
--databases '[path to fusion-report databases]'
```

### `--ericscript_ref`

Required reference in order to run `EricScript`.

```bash
--ericscript_ref '[path to EricScript reference]'
```

### `--fasta`

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--fasta '[path to Fasta reference]'
```

### `--fusioncatcher_ref`

Required reference in order to run `Fusioncatcher`.

```bash
--fusioncatcher_ref '[path to Fusioncatcher reference]'
```

### `--genome`

This pipeline uses only `Homo Sapiens` version `GRCh38`. Also make sure to specify `--genomes_base`.

```bash
--genome 'GRCh38' --genome_base '/path/to/references'
```

### `--gtf`

Required annotation file.

```bash
--gtf '[path to GTF annotation]'
```

### `--reference_release`

Ensembl version.

```bash
# ftp://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/
--reference_release '97'
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

### `--transcript`

Required transcript file.

```bash
--transcript '[path to transcript reference]'
```

## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack).

## AWS Batch specific parameters

Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use [`-profile awsbatch`](https://github.com/nf-core/configs/blob/master/conf/awsbatch.config) and then specify all of the following parameters.

### `--awsqueue`

The JobQueue that you intend to use on AWS Batch.

### `--awsregion`

The AWS region in which to run your job. Default is set to `eu-west-1` but can be adjusted to your needs.

### `--awscli`

The [AWS CLI](https://www.nextflow.io/docs/latest/awscloud.html#aws-cli-installation) path in your custom AMI. Default: `/home/ec2-user/miniconda/bin/aws`.

The AWS region to run your job in. Default is set to `eu-west-1` but can be adjusted to your needs.
Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

### `--debug`

To run only a specific tool (testing freshly implemented tool) just add `--debug` parameter. This parameter only works on **fusion tools only**!

### `--read_length`

Length is used to build a STAR index. Default is 100bp (Illumina).

### `--outdir`

The output directory where the results will be saved.

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `--email_on_fail`

This works exactly as with `--email`, except emails are only sent if the workflow is not successful.

### `--max_multiqc_email_size`

Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB).

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.
You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--custom_config_version`

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default: `master`.

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
