# nf-core/rnafusion: Usage <!-- omit in toc -->

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/rnafusion/usage](https://nf-co.re/rnafusion/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Pipeline summary

The pipeline is divided into two parts:

1. Download and build references

- specified with `--build_references` parameter
- required only once before running the pipeline
- **Important**: has to be run with each new release

2. Detecting fusions

- Supported tools: `Arriba`, `FusionCatcher`, `STAR-Fusion`, and `StringTie`
- QC: `Fastqc`, `MultiQC`, and `Picard CollectInsertSize`, `Picard CollectWgsMetrics`, `Picard Markduplicates`
- Fusions visualization: `Arriba`, `fusion-report`, `FusionInspector`, and `vcf_collect`

## Download and build references

The rnafusion pipeline needs references for the fusion detection tools, so downloading these is a **requirement**.

> **IMPORTANT**
>
> - Note that this step takes about 24 hours to complete on HPC.
> - Do not provide a samplesheet via the `input` parameter, otherwise the pipeline will run the analysis directly after downloading the references (except if that is what you want).

```bash
nextflow run nf-core/rnafusion \
  -profile <docker/singularity/.../institute> \
  --build_references --all \
  --cosmic_username <EMAIL> --cosmic_passwd <PASSWORD> \
  --genomes_base <PATH/TO/REFERENCES> \
  --outdir <PATH/TO/REFERENCES>
```

References for each tools can also be downloaded separately with:

```bash
nextflow run nf-core/rnafusion \
  -profile <docker/singularity/.../institute> \
  --build_references --<tool1> --<tool2> ... \
  --cosmic_username <EMAIL> --cosmic_passwd <PASSWORD> \
  --genomes_base <PATH/TO/REFERENCES> \
  --outdir <OUTPUT/PATH>
```

### Downloading the cosmic database with SANGER or QUIAGEN

#### For academic users

First register for a free account at COSMIC at [https://cancer.sanger.ac.uk/cosmic/register](https://cancer.sanger.ac.uk/cosmic/register) using a university email. The account is **only activated upon** clicking the link in the registration email.

#### For non-academic users

Use credentials from QIAGEN and add `--qiagen`

```bash
nextflow run nf-core/rnafusion \
  -profile <docker/singularity/.../institute> \
  --build_references --<tool1> --<tool2> ... \
  --cosmic_username <EMAIL> --cosmic_passwd <PASSWORD> \
  --genomes_base <PATH/TO/REFERENCES> \
  --outdir <OUTPUT/PATH> --qiagen
```

#### STAR-Fusion references downloaded vs built

By default STAR-Fusion references are **built**. You can also download them from [CTAT](https://github.com/NCIP/Trinity_CTAT/wiki) by using the flag `--starfusion_build FALSE` for both reference building and fusion detection. This allows more flexibility for different organisms but **be aware that STAR-Fusion reference download is not recommended as not fully tested!**

#### Issues with building references

If process `FUSIONREPORT_DOWNLOAD` times out, it could be due to network restriction (for example if trying to run on HPC). As this process is lightweight in cpu, memory and time, running on local machines with the following options might solve the issue:

```bash
nextflow run nf-core/rnafusion  \
  -profile <docker/singularity/.../institute> \
  --build_references \
  --cosmic_username <EMAIL> --cosmic_passwd <PASSWORD> \
  --fusionreport \
  --genomes_base <PATH/TO/REFERENCES> \
  --outdir <OUTPUT/PATH>
```

Adjustments for cpu and memory requirements can be done by feeding a custom configuration with `-c /PATH/TO/CUSTOM/CONFIG`.
Where the custom configuration could look like (adaptation to local machine necessary):

```text
process {
  withName:  'NFCORE_RNAFUSION:BUILD_REFERENCES:FUSIONREPORT_DOWNLOAD' {
    memory = '8.GB'
    cpus = 4
  }
}
```

The four `fusion-report` files: `cosmic.db`, `fusiongdb2.db`, `mitelman.db`
should then be copied into the HPC `<REFERENCE_PATH>/references/fusion_report_db`.

#### Note about fusioncatcher references

The references are only built based on ensembl version 102. It is not possible currently to use any other version/source.

## Running the pipeline

### Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. The pipeline will detect whether a sample is single- or paired-end from the samplesheet - the `fastq_2` column is empty for single-end. The samplesheet has to be a comma-separated file (.csv) but can have as many columns as you desire. There is a strict requirement for the first 4 columns to match those defined in the table below with the header row included.
A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice.

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,strandedness
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,forward
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz,forward
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz,forward
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,,forward
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,,forward
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,,forward
TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz,,forward
```

As you can see above for multiple runs of the same sample, the `sample` name has to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis.

| Column         | Description                                                                                                                                                                            |
| -------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`       | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1`      | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2`      | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `strandedness` | Strandedness: forward or reverse.                                                                                                                                                      |

### Starting commands

The pipeline can either be run using all fusion detection tools or specifying individual tools. Visualisation tools will be run on all fusions detected. To run all tools (`arriba`, `fusioncatcher`, `starfusion`, `stringtie`) use the `--all` parameter:

```bash
nextflow run nf-core/rnafusion \
  -profile <docker/singularity/.../institute> \
  --all \
  --input <SAMPLE_SHEET.CSV> \
  --genomes_base <PATH/TO/REFERENCES> \
  --outdir <OUTPUT/PATH>
```

To run only a specific detection tool use: `--tool`:

```bash
nextflow run nf-core/rnafusion \
  -profile <docker/singularity/.../institute> \
  --<tool1> --<tool2> ... \
  --input <SAMPLE_SHEET.CSV> \
  --genomes_base <PATH/TO/REFERENCES> \
  --outdir <OUTPUT/PATH>
```

> **IMPORTANT: Either `--all` or `--<tool>`** is necessary to run detection tools

`--genomes_base` should be the path to the directory containing the folder `references/` that was built with `--build_references`.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/rnafusion -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

:::warning
Conda is not currently supported.
Supported genome is currently only GRCh38.
:::

### Options

#### Trimming

When the flag `--fastp_trim` is used, `fastp` is used to provide all tools with trimmed reads. Quality and adapter trimming by default. In addition, tail trimming and adapter_fastq specification are possible. Example usage:

```bash
nextflow run nf-core/rnafusion \
-profile <docker/singularity/.../institute> \
--<tool1> --<tool2> ... \
--input <SAMPLE_SHEET.CSV> \
--genomes_base <PATH/TO/REFERENCES> \
--outdir <OUTPUT/PATH> \
--fastp_trim \
--trim_tail <INTEGER> (optional) \
--adapter_fastq <PATH/TO/ADAPTER/FASTQ> (optional)
```

#### Filter fusions detected by 2 or more tools

```bash
nextflow run nf-core/rnafusion \
  -profile <docker/singularity/.../institute> \
  --<tool1> --<tool2> ... \
  --input <SAMPLE_SHEET.CSV> \
  --genomes_base <PATH/TO/REFERENCES> \
  --outdir <OUTPUT/PATH>
  --tools_cutoff <INT>
```

`--tools_cutoff INT` will discard fusions detected by less than INT tools both for display in fusionreport html index and to consider in fusioninspector. Default = 1, no filtering.

#### Adding custom fusions to consider as well as the detected set: whitelist

```bash
nextflow run nf-core/rnafusion \
  -profile <docker/singularity/.../institute> \
  --<tool1> --<tool2> ... \
  --input <SAMPLE_SHEET.CSV> \
  --genomes_base <PATH/TO/REFERENCES> \
  --outdir <OUTPUT/PATH>
  --whitelist <WHITELIST/PATH>
```

The custom fusion file should have the following format:

```
GENE1--GENE2
GENE3--GENE4
```

#### Running FusionInspector only

FusionInspector can be run as a standalone with:

```bash
nextflow run nf-core/rnafusion \
-profile <docker/singularity/.../institute> \
--fusioninspector_only \
--fusioninspector_fusions <PATH_TO_CUSTOM_FUSION_FILE> \
--input <SAMPLE_SHEET.CSV> \
--outdir <PATH>
```

The custom fusion file should have the following format:

```
GENE1--GENE2
GENE3--GENE4
```

#### Skipping QC

```bash
nextflow run nf-core/rnafusion \
-profile <docker/singularity/.../institute> \
--skip_qc \
--all OR <--tool>
--input <SAMPLE_SHEET.CSV> \
--genomes_base <PATH/TO/REFERENCES> \
--outdir <PATH>
```

This will skip all QC-related processes (picard metrics collection)

#### Skipping visualisation

```bash
nextflow run nf-core/rnafusion \
-profile <docker/singularity/.../institute> \
--skip_vis \
--all OR <--tool>
--input <SAMPLE_SHEET.CSV> \
--genomes_base <PATH/TO/REFERENCES> \
--outdir <PATH>
```

This will skip all visualisation processes, including `fusion-report`, `FusionInspector` and `Arriba` visualisation.

#### Optional manual feed-in of fusion files

It is possible to give the output of each tool manually using the argument: `--<tool>_fusions PATH/TO/FUSION/FILE`: this feature need more testing, don't hesitate to open an issue if you encounter problems.

#### Set different `--limitSjdbInsertNsj` parameter

There are two parameters to increase the `--limitSjdbInsertNsj` parameter if necessary:

- `--fusioncatcher_limitSjdbInsertNsj`, default: 2000000
- `--fusioninspector_limitSjdbInsertNsj`, default: 1000000

Use the parameter `--cram` to compress the BAM files to CRAM for specific tools. Options: arriba, starfusion. Leave no space between options:

- `--cram arriba,starfusion`, default: []
- `--cram arriba`

### Troubleshooting

#### GstrandBit issues

The issue below sometimes occurs:

```
EXITING because of FATAL ERROR: cannot insert sequence on the fly because of strand GstrandBit problem
SOLUTION: please contact STAR author at https://groups.google.com/forum/#!forum/rna-star
```

As the error message suggests, it is a STAR-related error and your best luck in solving it will be the forum.

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/rnafusion
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/rnafusion releases page](https://github.com/nf-core/rnafusion/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
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
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.
- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
  - Needs to run in two steps: with `--build_references` first and then without `--build_references` to run the analysis
  - !!!! Run with `-stub` as all references need to be downloaded otherwise !!!!

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs). -->

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
