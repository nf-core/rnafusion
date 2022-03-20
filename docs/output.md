# nf-core/rnafusion: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview <!-- omit in toc -->

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* [Download and build references](#references) - Build references needed to run the rest of the pipeline
* [STAR](#star) - Alignment for arriba, squid and STAR-fusion
* [cat](#cat) - Concatenated fastq files per sample ID
* [Arriba](#arriba) - Arriba fusion detection
* [Pizzly](#pizzly) - Pizzly fusion detection
* [Squid](#squid) - Squid fusion detection
* [STAR-fusion](#starfusion) - STAR-fusion fusion detection
* [FusionCatcher](#fusioncatcher) - Fusion catcher fusion detection
* [Arriba visualisation](#arriba-visualisation) - Arriba visualisation report
* [Fusion-report](#fusion-report) - Summary of the findings of each tool and comparison to COSMIC, Mitelman and FusionGBD databases
* [FusionInspectior](#fusionInspector) - IGV-based visualisation tool for fusions filtered by fusion-report
* [FastQC](#fastqc) - Raw read QC
* [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
* [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution



<!-- 1. Input samplesheet check
2. Concatenate fastq files per sample
3. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))
5. Arriba subworkflow
    * STAR alignment
    * Samtool sort
    * Samtool index
    * [Arriba](https://github.com/suhrig/arriba) fusion detection
    * [Arriba](https://github.com/suhrig/arriba) visualisation
6. Pizzly subworkflow
    * [Kallisto](https://pachterlab.github.io/kallisto/) quantification
    * [Pizzly](https://github.com/pmelsted/pizzly) fusion detection
7. Squid subworkflow
    * [STAR](https://github.com/alexdobin/STAR) alignment
    * [Samtools view](http://www.htslib.org/): convert sam output from STAR to bam
    * [Samtools sort](http://www.htslib.org/): bam output from STAR
    * [Squid](https://github.com/Kingsford-Group/squid) fusion detection
    * [Squid](https://github.com/Kingsford-Group/squid) annotate
8. STAR-fusion subworkflow
    * [STAR](https://github.com/alexdobin/STAR) alignment
    * [STAR-fusion](https://github.com/STAR-Fusion/STAR-Fusion) fusion detection
9. Fusioncatcher subworkflow
    * [FusionCatcher](https://github.com/ndaniel/fusioncatcher) fusion detection
10. Fusion-report subworkflow
    * Merge all fusions detected by the different tools
    * [Fusion-report](https://github.com/matq007/fusion-report)
11. FusionInspector subworkflow
    * [FusionInspector](https://github.com/FusionInspector/FusionInspector) -->

<!-- TODO details on each step -->
### Download and build references

<details markdown="1">
<summary>Output files</summary>

 * `genomes_base/`
    * `arriba`

<!-- TODO complete list of reference output files -->


### STAR

<details markdown="1">
<summary>Output files</summary>

<!-- TODO complete list of output files -->

STAR is run 3 times:

For arriba with the parameters:

>`--readFilesCommand zcat \
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
--chimMultimapNmax 50`

For squid with the parameters:
>`--twopassMode Basic \
--chimOutType SeparateSAMold \
--chimSegmentMin 20 \
--chimJunctionOverhangMin 12 \
--alignSJDBoverhangMin 10 \
--outReadsUnmapped Fastx \
--outSAMstrandField intronMotif \
--outSAMtype BAM SortedByCoordinate \
--readFilesCommand zcat`

For STAR-fusion with the parameters:
>`--twopassMode Basic \
--outReadsUnmapped None \
--readFilesCommand zcat \
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
--chimOutType Junctions`



Alignment for arriba, squid and STAR-fusion
* [cat](#cat) - Concatenated fastq files per sample ID
* [Arriba](#arriba) - Arriba fusion detection
* [Pizzly](#pizzly) - Pizzly fusion detection
* [Squid](#squid) - Squid fusion detection
* [STAR-fusion](#starfusion) - STAR-fusion fusion detection
* [FusionCatcher](#fusioncatcher) - Fusion catcher fusion detection
* [Arriba visualisation](#arriba-visualisation) - Arriba visualisation report
* [Fusion-report](#fusion-report) - Summary of the findings of each tool and comparison to COSMIC, Mitelman and FusionGBD databases
* [FusionInspectior](#fusionInspector) - IGV-based visualisation tool for fusions filtered by fusion-report

### FastQC

<details markdown="1">
<summary>Output files</summary>

* `fastqc/`
    * `*_fastqc.html`: FastQC report containing quality metrics.
    * `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

* `multiqc/`
    * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
    * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
    * `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

* `pipeline_info/`
    * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
    * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
    * Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
