# nf-core/rnafusion: Usage <!-- omit in toc -->

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
  - [-profile](#-profile)
  - [-resume](#-resume)
  - [-c](#-c)
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

## Running the pipeline

The typical command for running the pipeline is as follows.

```bash
nextflow run nf-core/rnafusion \
  -profile docker \
  --input "*_R{1,2}.fastq.gz" \
  --arriba \
  --star_fusion \
  --fusioncatcher \
  --ericscript \
  --pizzly \
  --squid \
  --arriba_vis \
  --fusion_inspector
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

First, go to the [nf-core/rnafusion releases page](https://github.com/nf-core/rnafusion/releases) and find the latest version number - numeric only (eg. `1.2`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.2`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### -profile

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time.
For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

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

### -resume

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### -c

Specify the path to a specific config file (this is a core Nextflow command).
See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

#### Custom resource requests

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

Whilst these default requirements will hopefully work for most people with most data, you may find that you want to customise the compute resources that the pipeline requests. You can do this by creating a custom config file. For example, to give the workflow process `VEP` 32GB of memory, you could use the following config:

```nextflow
process {
  withName: VEP {
    memory = 32.GB
  }
}
```

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

### Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

#### Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
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

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.