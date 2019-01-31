# nf-core/rnafusion: Usage

## Table of contents

* [Introduction](#general-nextflow-info)
* [Running the pipeline](#running-the-pipeline)
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
* [Tool flags](#tool-flags)
  * [`--star_fusion`](#--star_fusion)
  * [`--fusioncatcher`](#--fusioncatcher)
    * [`--fc_extra_options`](#--fc_extra_options)
  * [`--ericscript`](#--ericscript)
  * [`--pizzly`](#--pizzly)  
  * [`--squid`](#--squid)
  * [`--tool_cutoff`](#--tool-cutoff)
  * [`--test`](#--test)
* [Reference Genomes](#reference-genomes)
  * [`--genome`](#--genome)
  * [`--fasta`](#--fasta)
  * [`--gtf`](#--gtf)
  * [`--star_index`](#--star_index)
  * [`--star_fusion_ref`](#--star_fusion_ref)
  * [`--fusioncatcher_ref`](#--fusioncatcher_ref)
  * [`--ericscript_ref`](#--ericscript_ref)
  * [`--pizzly_fasta`](#--pizzly_fasta)
  * [`--pizzly_gtf`](#--pizzly_gtf)
* [Job Resources](#job-resources)
* [Automatic resubmission](#automatic-resubmission)
* [Custom resource requests](#custom-resource-requests)
* [AWS batch specific parameters](#aws-batch-specific-parameters)
  * [`-awsbatch`](#-awsbatch)
    * [`--awsqueue`](#--awsqueue)
    * [`--awsregion`](#--awsregion)
* [Other command line parameters](#other-command-line-parameters)
  * [`--read_length`](#--read_length)
  * [`--singleEnd`](#--singleend)
  * [`--outdir`](#--outdir)
  * [`--email`](#--email)
  * [`-name`](#-name-single-dash)
  * [`-resume`](#-resume-single-dash)
  * [`-c`](#-c-single-dash)
  * [`--max_memory`](#--max_memory)
  * [`--max_time`](#--max_time)
  * [`--max_cpus`](#--max_cpus)
  * [`--plaintext_emails`](#--plaintext_emails)
  * [`--sampleLevel`](#--sampleLevel)
  * [`--multiqc_config`](#--multiqc_config)

## General Nextflow info

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline

The typical command for running the whole pipeline is as follows:

```bash
nextflow run nf-core/rnafusion --reads '*_R{1,2}.fastq.gz' --genome GRCh38 -profile docker --star_fusion --fusioncatcher --ericscript --pizzly --squid --fusion_inspector
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Is is also possible to execute specific tools:

```bash
nextflow run nf-core/rnafusion --reads '*_R{1,2}.fastq.gz' --genome GRCh38 -profile docker --fusioncatcher --ericscript
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

## Main Arguments

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile docker` - the order of arguments is important!

* `awsbatch`
  * A generic configuration profile to be used with AWS Batch.
* `conda`
  * A generic configuration profile to be used with [conda](https://conda.io/docs/)
* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
    * Pulls software from dockerhub: [`nfcore/rnafusion`](http://hub.docker.com/r/nfcore/rnafusion/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
    * Pulls software from docker-hub
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

## Tool flags

### `--star_fusion`

If enabled, executes `STAR-Fusion` tool.

### `--fusioncatcher`

If enabled, executes `Fusioncatcher` tool.

* `--fc_extra_options`
  * Parameter for specifying additional parameters. For more info, please refer to the documentation of the tool.

### `--ericscript`

If enabled, executes `Ericscript` tool.

### `--pizzly`

If enabled, executes `Pizzly` tool.

### `--squid`

If enabled, executes `Squid` tool.

### `--tool_cutoff`

This is a summary-report parameter which serves as a filter for how many tools are required to detect a fusion. The default number is 2. If a fusion is not detected by at least 2 tools it will not be displayed in the final report. Secondly, if number of tools used in the pipeline is less than 2, this filter is being ignored.

### `--test`

To run only a specific tool (testing freshly implemented tool) just add `--test` parameter. This parameter only works on **fusion tools only**!

```bash
nextflow run nf-core/rnafusion --reads '*_R{1,2}.fastq.gz' --genome GRCh38 -profile docker --star_fusion --test
```

## Reference Genomes

The pipeline config files come bundled with paths to the illumina iGenomes reference index files. If running with docker or AWS, the configuration is set up to use the [AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/) resource.

### `--genome` (using iGenomes)

There are 31 different species supported in the iGenomes references. To run the pipeline, you must specify which to use with the `--genome` flag.
You can find the keys to specify the genomes in the [iGenomes config file](../conf/igenomes.config). Common genomes that are supported are:

* Human
  * `--genome GRCh38` (recommended)

> **TL;DR** For now the pipeline only supports Human genome GRCh38.

Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the iGenomes resource. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to save such a file.

The syntax for this reference configuration is as follows:

```nextflow
params {
  genomes {
    'GRCh37' {
      fasta   = '<path to the genome fasta file>' // Used if no star index given
    }
    // Any number of additional genomes, key is used with --genome
  }
}
```

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

## Job Resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files in [`conf`](../conf) for examples.

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

### `--singleEnd`

By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--singleEnd` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--singleEnd --reads '*.fastq.gz'
```

### `--outdir`

The output directory where the results will be saved.

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to speicfy this on the command line for every run.

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

Note - you can use this to override defaults. For example, you can specify a config file using `-c` that contains the following:

```nextflow
process.$multiqc.module = []
```

### `--max_memory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'``

### `--max_time`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`

Set to receive plain-text e-mails instead of HTML formatted.

### `--multiqc_config`

Specify a path to a custom MultiQC configuration file.
