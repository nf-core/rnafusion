# nf-core/rnafusion: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview <!-- omit in toc -->

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Download and build references](#references) - Build references needed to run the rest of the pipeline
- [STAR](#star) - Alignment for arriba, squid and STAR-fusion
- [Cat](#cat) - Concatenated fastq files per sample ID
- [Arriba](#arriba) - Arriba fusion detection
- [Pizzly](#pizzly) - Pizzly fusion detection
- [Squid](#squid) - Squid fusion detection
- [STAR-fusion](#starfusion) - STAR-fusion fusion detection
- [FusionCatcher](#fusioncatcher) - Fusion catcher fusion detection
- [Samtools](#samtools) - SAM/BAM file manipulation
- [Arriba visualisation](#arriba-visualisation) - Arriba visualisation report
- [Fusion-report](#fusion-report) - Summary of the findings of each tool and comparison to COSMIC, Mitelman and FusionGBD databases
- [FusionInspector](#fusionInspector) - IGV-based visualisation tool for fusions filtered by fusion-report
- [Qualimap](#qualimap) - Quality control of alignment
- [Picard](#picard) - Collect metrics
- [FastQC](#fastqc) - Raw read quality control
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Download and build references

<details markdown="1">
<summary>Output files</summary>

- `genomes_base/`
  - `arriba`
    - `blacklist_hg38_GRCh38_v2.1.0.tsv.gz`
    - `protein_domains_hg38_GRCh38_v2.1.0.gff3`
    - `cytobands_hg38_GRCh38_v2.1.0.tsv`
  - `ensembl`
    - `Homo_sapiens.GRCh38.{ensembl_version}.all.fa`
    - `Homo_sapiens.GRCh38.{ensembl_version}.cdna.all.fa.gz`
    - `Homo_sapiens.GRCh38.{ensembl_version}.gtf`
    - `Homo_sapiens.GRCh38.{ensembl_version}.chr.gtf`
    - `Homo_sapiens.GRCh38.{ensembl_version}.chr.gtf.refflat`
  - `fusioncatcher`
    - `human_v<version>` - dir with all references for fusioncatcher
  - `fusion_report_db`
    - `cosmic.db`
    - `fusiongdb.db`
    - `fusiongdb2.db`
    - `mitelman.db`
  - `pizzly`
    - `kallisto` - file containing the kallisto index
  - `star` - dir with STAR index
  - `starfusion`
    - files and dirs used to build the index
    - `ctat_genome_lib_build_dir` - dir containing the index

(Only files or folders used by the pipeline are mentioned explicitly.)

</details>

### STAR

STAR is used to align to genome reference

STAR is run 3 times:

For arriba with the parameters:

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

For squid with the parameters:

```bash
--twopassMode Basic \
--chimOutType SeparateSAMold \
--chimSegmentMin 20 \
--chimJunctionOverhangMin 12 \
--alignSJDBoverhangMin 10 \
--outReadsUnmapped Fastx \
--outSAMstrandField intronMotif \
--outSAMtype BAM SortedByCoordinate \
--readFilesCommand zcat
```

For STAR-fusion with the parameters:

```bash
--twopassMode Basic \
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
--chimOutType Junctions
```

> STAR_FOR_STARFUSION uses `${params.ensembl_ref}/Homo_sapiens.GRCh38.${params.ensembl_version}.chr.gtf` whereas STAR_FOR_ARRIBA and STAR_FOR_SQUID use `${params.ensembl_ref}/Homo_sapiens.GRCh38.${params.ensembl_version}.gtf`

<details markdown="1">
<summary>Output files</summary>

- `star_for_<tool>`
_ **Common**
_ `<sample_id>.Log.final.out`
_ `<sample_id>.Log.progress.out`
_ `<sample_id>.SJ.out.tab`
_ **For arriba:**
_ `<sample_id>.Aligned.out.bam`
_ **For squid:**
_ `<sample_id>.Aligned.sortedByCoord.out.bam`
_ `<sample_id>.Chimeric.out.sam`
_ `<sample_id>.unmapped_1.fastq.gz`
_ `<sample_id>.unmapped_2.fastq.gz`
_ **For starfusion:**
_ `<sample_id>.Aligned.sortedByCoord.out.bam`
_ `<sample_id>.Chimeric.out.junction`
</details>

### Cat

Cat is used to concatenate fastq files belonging to the same sample.

<details markdown="1">
<summary>Output files</summary>

- `cat`
  - `<sample_id>_1.merged.fastq.gz`
  - `<sample_id>_2.merged.fastq.gz`

</details>

### Arriba

Arriba is used for i) detect fusion and ii) output a PDF report for the fusions found (visualisation):

#### Detection

<details markdown="1">
<summary>Output files</summary>

- `arriba`
  - `<sample_id>.arriba.fusions.tsv` - contains the identified fusions
  - `<sample_id>.arriba.fusions.discarded.tsv`

</details>

#### Visualisation

<details markdown="1">
<summary>Output files</summary>

- `arriba_visualisation`
  - `<sample_id>.pdf`

</details>

### Pizzly

The first step of the pizzly workflow is to run `kallisto quant`:

#### Kallisto

<details markdown="1">
<summary>Output files</summary>

- `kallisto`
  - `<sample_id>.kallisto_quant.fusions.txt`

</details>

Pizzly refines kallisto output.

#### Pizzly

Pizzly uses the following arguments:

```bash
-k 31 \
--align-score 2 \
--insert-size 400 \
--cache index.cache.txt
```

<details markdown="1">
<summary>Output files</summary>

- `pizzly`
  - `<sample_id>.pizzly.txt` - contains the identified fusions
  - `<sample_id>.pizzly.unfiltered.json`

</details>

### Squid

Squid is run in two steps: i) fusion detection and ii) fusion annotation but the output is in a common `squid` directory.

<details markdown="1">
<summary>Output files</summary>

- `squid`
  - `<sample_id>.squid.fusions_sv.txt` - contains the identified fusions
  - `<sample_id>.squid.fusions.annotated.txt`- contains the identified fusions annotatedvi

</details>

### STAR-fusion

<details markdown="1">
<summary>Output files</summary>

- `starfusion`
  - `<sample_id>.starfusion.fusion_predictions.tsv` - contains the identified fusions
  - `<sample_id>.starfusion.abridged.tsv`
  - `- contains the identified fusions.starfusion.abridged.coding_effect.tsv`

</details>

### FusionCatcher

<details markdown="1">
<summary>Output files</summary>

- `fusioncatcher`
_ `<sample_id>.fusioncatcher.fusion-genes.txt`
_ `<sample_id>.fusioncatcher.summary.txt` \* `<sample_id>.fusioncatcher.log`
</details>

### Samtools

#### Samtools view

Samtools view is used to convert the chimeric SAM output from STAR_FOR_SQUID to BAM

<details markdown="1">
<summary>Output files</summary>

- `samtools_view_for_squid`
  - `<sample_id>_chimeric.bam` - sorted BAM file

</details>

#### Samtools sort

Samtools sort is used to sort BAM files from STAR_FOR_ARRIBA (for arriba visualisation) and the chimeric BAM from STAR_FOR_SQUID

<details markdown="1">
<summary>Output files</summary>

- `samtools_sort_for_<arriba/squid>`
  - `<sample_id>(_chimeric)_sorted.bam` - sorted BAM file

</details>

#### Samtools index

Samtools index is used to index BAM files from STAR_FOR_ARRIBA (for arriba visualisation) and STAR_FOR_STARFUSION (for QC)

<details markdown="1">
<summary>Output files</summary>

- `samtools_for_<arriba/qc>`
  - `<sample_id>.(Aligned.sortedByCoord).out.bam.bai` -

</details>

### Fusion-report

<details markdown="1">
<summary>Output files</summary>

- `fusionreport`
  - <sample_id>
    - `<sample_id>.fusionreport.tsv`
    - `<sample_id>.fusionreport_filtered.tsv`
    - `index.html` - general report for all filtered fusions
    - `<fusion>.html` - specific report for each filtered fusion

</details>

### FusionInspector

<details markdown="1">
<summary>Output files</summary>

- `fusioninspector`
  - `<sample_id>.fusion_inspector_web.html` - visualisation report described in details [here](https://github.com/FusionInspector/FusionInspector/wiki/FusionInspector-Visualizations)
  - `FusionInspector.log`
  - `<sample_id>.FusionInspector.fusions.abridged.tsv`

</details>

### Qualimap

<details markdown="1">
<summary>Output files</summary>

- `qualimap`
  - `qualimapReport.html` - HTML report
  - `rnaseq_qc_results.txt` - TXT results
  - `css` - dir for html style
  - `images_qualimapReport`- dir for html images
  - `raw_data_qualimapReport` - dir for html raw data

</details>

### Picard

Picard CollectRnaMetrics and picard MarkDuplicates share the same outpur directory.

<details markdown="1">
<summary>Output files</summary>

- `picard`
  - `<sample_id>.MarkDuplicates.metrics.txt` - metrics from CollectRnaMetrics
  - `<sample_id>_rna_metrics.txt` - metrics from MarkDuplicates
  - `<sample_id>.bam` - BAM file with marked duplicates

</details>

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
