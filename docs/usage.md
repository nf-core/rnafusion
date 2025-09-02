# nf-core/rnafusion: Usage <!-- omit in toc -->

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/rnafusion/usage](https://nf-co.re/rnafusion/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Pipeline summary

The pipeline is divided into two parts:

1. Download and build references

- specified with `--references_only` parameter
- required only once before running the pipeline
- **Important**: has to be run with each new release

2. Detecting fusions

- Supported tools: `Arriba`, `FusionCatcher`, `STAR-Fusion`, `StringTie` and `CTAT-SPLICING`
- QC: `Fastqc`, `MultiQC`, and `Picard CollectInsertSize`, `Picard CollectWgsMetrics`, `Picard Markduplicates`
- Fusions visualization: `Arriba`, `fusion-report`, `FusionInspector`, and `vcf_collect`

## Download and build references

The rnafusion pipeline needs references for the fusion detection tools, so downloading these is a **requirement**.

The references for the pipeline can be downloaded from the nf-core AWS megatests S3 bucket using the following command for the [AWS CLI tool](https://github.com/aws/aws-cli):

```bash
aws --no-sign-request s3 sync s3://nf-core-awsmegatests/rnafusion/references/ <path_to_references>
```

The path to the downloaded references can then be provided to the pipeline with the `--genomes_base` parameter.

:warning: **Please note that the references are large and can take a long time to download, so it is recommended to download them once and use them for all future runs of the pipeline.**

The fusion report references available on the S3 bucket do not contain information from cosmic due to licensing issues. If you want to use the cosmic database, you will need to build the fusion report references yourself by either deleting the `fusion_report_db` folder in the references directory or by specifying a different location for the fusion report directory with `--fusionreport_ref <PATH/TO/REFERENCES>`. The cosmic username and password should also be given in this case using `--cosmic_username <EMAIL>` and `--cosmic_passwd <PASSWORD>`.

Additionally, the references can be built by the pipeline using the following command:

```bash
nextflow run nf-core/rnafusion \
  -profile <docker/singularity/.../institute> \
  --references_only --tools all \
  --cosmic_username <EMAIL> --cosmic_passwd <PASSWORD> \
  --genomes_base <PATH/TO/REFERENCES> \
  --outdir <PATH/TO/REFERENCES>
```

> **IMPORTANT**
>
> - Note that this step takes about 24 hours to complete on HPC.

References for each tools can also be downloaded separately with:

```bash
nextflow run nf-core/rnafusion \
  -profile <docker/singularity/.../institute> \
  --references_only --tools <comma-separated-list-of-tools> \
  --cosmic_username <EMAIL> --cosmic_passwd <PASSWORD> \
  --genomes_base <PATH/TO/REFERENCES> \
  --outdir <OUTPUT/PATH>
```

If you are not covered by the research COSMIC license and want to avoid using COSMIC, you can provide the additional option `--no_cosmic`.

### Downloading the cosmic database with SANGER or QIAGEN

#### For academic users

First register for a free account at COSMIC at [https://cancer.sanger.ac.uk/cosmic/register](https://cancer.sanger.ac.uk/cosmic/register) using a university email. The account is **only activated upon** clicking the link in the registration email.

#### For non-academic users

Use credentials from QIAGEN and add `--qiagen`

```bash
nextflow run nf-core/rnafusion \
  -profile <docker/singularity/.../institute> \
  --references_only --tools <comma-separated-list-of-tools> \
  --cosmic_username <EMAIL> --cosmic_passwd <PASSWORD> \
  --genomes_base <PATH/TO/REFERENCES> \
  --outdir <OUTPUT/PATH> --qiagen
```

<!-- #### STAR-Fusion references downloaded vs built

By default STAR-Fusion references are **built**. You can also download them from [CTAT](https://github.com/NCIP/Trinity_CTAT/wiki) by using the flag `--starfusion_build FALSE` for both reference building and fusion detection. This allows more flexibility for different organisms but **be aware that STAR-Fusion reference download is not recommended as not fully tested!** -->

#### Issues with building references

If process `FUSIONREPORT_DOWNLOAD` times out, it could be due to network restriction (for example if trying to run on HPC). As this process is lightweight in cpu, memory and time, running on local machines with the following options might solve the issue:

```bash
nextflow run nf-core/rnafusion  \
  -profile <docker/singularity/.../institute> \
  --references_only \
  --cosmic_username <EMAIL> --cosmic_passwd <PASSWORD> \
  --tools fusionreport \
  --genomes_base <PATH/TO/REFERENCES> \
  --outdir <OUTPUT/PATH>
```

Adjustments for cpu and memory requirements can be done by feeding a custom configuration with `-c /PATH/TO/CUSTOM/CONFIG`.
Where the custom configuration could look like (adaptation to local machine necessary):

```text
process {
  withName:  'NFCORE_RNAFUSION:RNAFUSION:BUILD_REFERENCES:FUSIONREPORT_DOWNLOAD' {
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

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. The pipeline will detect whether a sample is single- or paired-end from the samplesheet - the `fastq_2` column is empty for single-end. The samplesheet has to be a comma-separated (.csv), tab-separated (.tsv), yaml (.yaml or .yml) or json (.json) file but can have as many columns as you desire. There is a strict requirement for the `sample` and `strandedness` columns. One or more of these columns should be provided too: `fastq_1`, `bam`, `cram`, `junctions` and `splice_junctions`
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

As you can see above for multiple runs of the same sample, the `sample` name has to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Note that multiple rows per sample are not supported for samples that contain `bam`, `cram`, `junctions` and/or `splice_junctions` files.

| Column             | Description                                                                                                                                                                                                                                                                                                                                                                                                                                   | Required           |
| ------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------ |
| `sample`           | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`).                                                                                                                                                                                                                                                        | :white_check_mark: |
| `strandedness`     | Strandedness: forward or reverse.                                                                                                                                                                                                                                                                                                                                                                                                             | :white_check_mark: |
| `fastq_1`          | Full path to FastQ file for Illumina short reads 1. File must exist, has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". It's recommended to always provide the FASTQ files because the pipeline will be able to create any missing files from these. The FASTQ files are required to run `salmon`, `fusioninspector` and `fusioncatcher`.                                                                                      | :grey_question:    |
| `fastq_2`          | Full path to FastQ file for Illumina short reads 2. File must exist, has to be gzipped and have the extension ".fastq.gz" or ".fq. It's recommended to always provide the FASTQ files because the pipeline will be able to create any missing files from these. The FASTQ files are required to run `salmon`, `fusioninspector` and `fusioncatcher`.gz".                                                                                      | :grey_question:    |
| `bam`              | Full path to the BAM file created with STAR. File has to exist and must have the extension ".bam". It's the responsibility of the pipeline user to make sure this file has been correctly created, see the [prepare chapter](#preparing-bamcramjunctionssplice_junctions) for more information. The BAM file is required to run `ctatsplicing`, `stringtie`, `fusioninspector` and `arriba` when the `fastq_1` and `cram` fields are empty.   | :grey_question:    |
| `bai`              | Full path to the index of the BAM file. File has to exist and must have the extension ".bai".                                                                                                                                                                                                                                                                                                                                                 | :x:                |
| `cram`             | Full path to the CRAM file created with STAR. File has to exist and must have the extension ".cram". It's the responsibility of the pipeline user to make sure this file has been correctly created, see the [prepare chapter](#preparing-bamcramjunctionssplice_junctions) for more information. The CRAM file is required to run `ctatsplicing`, `stringtie`, `fusioninspector` and `arriba` when the `fastq_1` and `bam` fields are empty. | :grey_question:    |
| `crai`             | Full path to the index of the CRAM file. File has to exist and must have the extension ".crai".                                                                                                                                                                                                                                                                                                                                               | :x:                |
| `junctions`        | Full path to the file containing chimeric junctions determined by STAR. File has to exist and must have the extension ".junction". It's the responsibility of the pipeline user to make sure this file has been correctly created, see the [prepare chapter](#preparing-bamcramjunctionssplice_junctions) for more information. The junctions file is required to run `starfusion` and `ctatsplicing` when the `fastq_1` field is empty.      | :grey_question:    |
| `splice_junctions` | Full path to the file containing splice junctions determined by STAR. File has to exist and must have the extension ".SJ.out.tab". It's the responsibility of the pipeline user to make sure this file has been correctly created, see the [prepare chapter](#preparing-bamcramjunctionssplice_junctions) for more information. The splice junctions file is required to run `ctatsplicing` when the `fastq_1` field is empty.                | :grey_question:    |
| `seq_platform`     | The sequencing platform used to create to sequence the data in the FASTQ files. This value will take precedence over the value provided with `--seq_platform`.                                                                                                                                                                                                                                                                                | :x:                |
| `seq_center`       | The sequencing center in which the data in the FASTQ files was sequenced. This value will take precedence over the value provided with `--seq_center`.                                                                                                                                                                                                                                                                                        | :x:                |

:white_check_mark: = Required
:x: = Not required
:grey_question: = One of these columns should be provided

### Preparing BAM/CRAM/junctions/splice_junctions

In the pipeline the following STAR command is used to produce the needed files:

```bash
STAR \\
    --genomeDir <path-to-star-index> \
    --readFilesIn <comma-separated-list-of-forward-fastqs> <comma-separated-list-of-reverse-fastqs> \
    --runThreadN <threads> \
    --outFileNamePrefix <sample-name>. \
    --outSAMattrRGline 'ID:<sample-name>' 'SM:<sample-name>' \
    --outReadsUnmapped None  \
    --outSAMstrandField intronMotif \
    --chimOutJunctionFormat 1 \
    --twopassMode None \
    --outFilterMultimapNmax 50 \
    --chimMultimapNmax 50 \
    --quantMode GeneCounts \
    --outSAMunmapped Within \
    --readFilesCommand zcat  \
    --alignSJstitchMismatchNmax 5 -1 5 5 \
    --outSAMtype BAM SortedByCoordinate \
    --chimSegmentMin 10 \
    --peOverlapNbasesMin 10 \
    --alignSplicedMateMapLminOverLmate 0.5 \
    --chimJunctionOverhangMin 10 \
    --chimScoreJunctionNonGTAG 0 \
    --chimScoreDropMax 30 \
    --chimScoreSeparation 1  \
    --chimSegmentReadGapMax 3 \
    --chimOutType Junctions WithinBAM'
```

We found that this command produces the best results for all downstream processes in the pipeline. It is highly recommended to use the same command for the input BAM, CRAM, junctions and splice_junctions files.

The pipeline will still work when another command has been used, but it is possible that the results will be significantly different from the standard flow.

### Starting commands

The pipeline can either be run using all fusion detection tools or by specifying individual tools. Visualisation tools will be run on all fusions detected. To run all tools (`arriba`, `ctatsplicing`, `fusioncatcher`, `starfusion`, `stringtie`, `fusionreport`, `fastp`, `salmon`, `fusioninspector`) use the `all` option for the `--tools` parameter:

```bash
nextflow run nf-core/rnafusion \
  -profile <docker/singularity/.../institute> \
  --tools all \
  --input <SAMPLE_SHEET.CSV> \
  --genomes_base <PATH/TO/REFERENCES> \
  --outdir <OUTPUT/PATH>
```

To run only a set of specific tools use `--tools` with a comma separated list of requested tools:

```bash
nextflow run nf-core/rnafusion \
  -profile <docker/singularity/.../institute> \
  --tools <comma-separated-list-of-tools> \
  --input <SAMPLE_SHEET.CSV> \
  --genomes_base <PATH/TO/REFERENCES> \
  --outdir <OUTPUT/PATH>
```

If you are not covered by the research COSMIC license and want to avoid using COSMIC, you can provide the additional option `--no_cosmic`.

> **IMPORTANT: `--tools`** is necessary to run detection tools

`--genomes_base` should be the path to the directory containing the folder `references/` that was built with `--references_only`.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

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

When the flag `fastp` tool is used in `--tools`, `fastp` is used to provide all tools with trimmed reads. Quality and adapter trimming by default. In addition, tail trimming and adapter_fastq specification are possible. Example usage:

```bash
nextflow run nf-core/rnafusion \
-profile <docker/singularity/.../institute> \
--tools fastp,... \
--input <SAMPLE_SHEET.CSV> \
--genomes_base <PATH/TO/REFERENCES> \
--outdir <OUTPUT/PATH> \
--fastp_trim \
--trim_tail <INTEGER> (optional) \
--adapter_fastq <PATH/TO/ADAPTER/FASTQ> (optional)
```

The additional `--trim_tail_fusioncatcher` flag will toggle an additional `fastp` process, especially useful is reads are above 100 bp, which is not handled well by FusionCatcher. The parameter `--trim_tail_fusioncatcher`need to be provided with the number of bases to remove from the tail end.

#### Filter fusions detected by 2 or more tools

```bash
nextflow run nf-core/rnafusion \
  -profile <docker/singularity/.../institute> \
  --tools <comma-separated-list-of-tools> \
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
  --tools <comma-separated-list-of-tools> \
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
--tools fusioninspector \
--skip_qc \
--skip_vis \
--skip_vcf \
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
--tools <comma-separated-list-of-tools> \
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
--tools <comma-separated-list-of-tools> \
--input <SAMPLE_SHEET.CSV> \
--genomes_base <PATH/TO/REFERENCES> \
--outdir <PATH>
```

This will skip all visualisation processes, including `fusion-report`, `FusionInspector` and `Arriba` visualisation.

#### Optional manual feed-in of fusion files

It is possible to give the output of each fusion detection tool manually using the argument: `--<tool>_fusions PATH/TO/FUSION/FILE`: this feature needs more testing, don't hesitate to open an issue if you encounter problems.

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

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/rnafusion releases page](https://github.com/nf-core/rnafusion/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

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
  - !!!! Run with `-stub` as all references need to be downloaded otherwise !!!!

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

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
