# nf-core/rnafusion: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview <!-- omit in toc -->

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Download and build references](#references) - Build references needed to run the rest of the pipeline
- [STAR](#star) - Alignment for arriba, and STAR-fusion
- [Cat](#cat) - Concatenate fastq files per sample ID
- [Arriba](#arriba) - Arriba fusion detection
- [STAR-fusion](#starfusion) - STAR-fusion fusion detection
- [StringTie](#stringtie) - StringTie assembly
- [FusionCatcher](#fusioncatcher) - Fusion catcher fusion detection
- [Samtools](#samtools) - SAM/BAM file manipulation
- [Fusion-report](#fusion-report) - Summary of the findings of each tool and comparison to COSMIC, Mitelman, and FusionGDB2 databases
- [FusionInspector](#fusionInspector) - Supervised analysis of fusion predictions from fusion-report, recover and re-score evidence for such predictions
- [Arriba visualisation](#arriba-visualisation) - Arriba visualisation report for FusionInspector fusions
- [Picard](#picard) - Collect QC metrics
- [FastQC](#fastqc) - Raw read quality control
- [Salmon](#salmon) - Normalized gene expression calculation
- [MultiQC](#multiqc) - Aggregate reports describing QC results from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## Download and build references

<details markdown="1">
<summary>Output reference files and folder structure</summary>

### References directory structure

- `references/`
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
    - `Homo_sapiens.GRCh38.{ensembl_version}.interval_list`
  - `fusioncatcher`
    - `human_v<version>` - dir with all references for fusioncatcher
  - `fusion_report_db`
    - `cosmic.db`
    - `fusiongdb2.db`
    - `mitelman.db`
  - `star` - dir with STAR index
  - `starfusion`
    - files and dirs used to build the index
    - `ctat_genome_lib_build_dir` - dir containing the index

(Only files or folders used by the pipeline are mentioned explicitly.)

</details>

## Main pipeline workflow

> If no argument is specified here, the tool was used with default parameters.

### Directory structure

```text
{outdir}
├── arriba
├── arriba_visualisation
├── cram_arriba
├── cram_starfusion
├── fastp
├── fastqc
├── fusioncatcher
├── fusioninspector
├── fusionreport
├── kallisto_quant
├── megafusion
├── multiqc
├── picard
├── pipeline_info
├── samtools_sort_for_arriba
├── star_for_arriba
├── star_for_starfusion
├── starfusion
└── work
.nextflow.log
```

If no parameters are specified, the default is applied.

### Arriba

[Arriba](https://arriba.readthedocs.io/en/latest/) is used for i) detect gene fusions and ii) create a PDF report for the fusions found (visualisation):

#### Detection

<details markdown="1">
<summary>Output files</summary>

- `arriba/`
  - `<sample>.arriba.fusions.tsv` - contains the identified fusions
  - `<sample>.arriba.fusions.discarded.tsv`

</details>

#### Visualisation

<details markdown="1">
<summary>Output files</summary>

- `arriba_visualisation/`
  - `<sample>_combined_fusions_arriba_visualisation.pdf`

</details>

The visualisation displays the fusions that fusioninspector outputs. That means that fusions from all callers are aggregated (by fusion-report) and then analyzed through fusioninspector (Note: Fusioninspecor contains a filtering step!).

### Cat

<details markdown="1">
<summary>Output files</summary>

- `cat/`
  - `<sample>_1.merged.fastq.gz`
  - `<sample>_2.merged.fastq.gz`

</details>

If multiple libraries or runs have been provided for the same sample in the input samplesheet (e.g. to increase sequencing depth) then these will be merged at the very beginning of the pipeline in order to have consistent sample naming throughout the pipeline. Please refer to the [usage](https://nf-co.re/rnafusion/usage#samplesheet-input) documentation to see how to specify these samples in the input samplesheet.

### Fastp

If `--trim_fastp` is selected, [fastp](https://github.com/OpenGene/fastp) will filter low quality reads as well as bases at the 5' and 3' ends, trim adapters (automatically detected, but input with parameter `--adapter_fasta` is possible). 3' trimming is also possible via parameter `--trim_tail`.

<details markdown="1">
<summary>Output files</summary>

- `fastp/`
  - `<sample>_1.fastp.fastq.gz`
  - `<sample>_2.fastp.fastq.gz`
  - `<sample>.fastp.html`
  - `<sample>.fastp.json`
  - `<sample>.fastp.log`

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

:::note
The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.
:::

### FusionCatcher

<details markdown="1">
<summary>Output files</summary>

- `fusioncatcher`
  - `<sample>.fusioncatcher.fusion-genes.txt`
  - `<sample>.fusioncatcher.summary.txt`
  - `<sample>.fusioncatcher.log`

</details>

[FusionCatcher](https://github.com/ndaniel/fusioncatcher) searches for novel/known somatic fusion genes translocations, and chimeras in RNA-seq data. Possibility to use parameter `--fusioncatcher_limitSjdbInsertNsj` to modify limitSjdbInsertNsj.

### FusionInspector

<details markdown="1">
<summary>Output files</summary>

- `fusioninspector`
  - `<sample>.fusion_inspector_web.html` - visualisation report described in details [here](https://github.com/FusionInspector/FusionInspector/wiki/FusionInspector-Visualizations)
  - `FusionInspector.log`
  - `<sample>.FusionInspector.fusions.abridged.tsv`

</details>

[FusionInspector](https://github.com/FusionInspector/FusionInspector/tree/master) performs a validation of fusion transcript predictions. Possibility to use `--fusioninspector_limitSjdbInsertNsj` to set limitSjdbInsertNsj to anything other than the default 1000000.

### Fusion-report

Please note that fusion-report is executed from fork https://github.com/Clinical-Genomics/fusion-report

<details markdown="1">
<summary>Output files</summary>

- `fusionreport`
  - <sample>
    - `<sample>.fusionreport.tsv`
    - `<sample>.fusionreport_filtered.tsv`
    - `<sample>_fusionreport_index.html` - general report for all filtered fusions
    - `<sample>.fusions.csv` - index in csv format
    - `<fusion>_<pair>.html` - specific report for each filtered fusion

</details>

[Fusion-report](https://github.com/matq007/fusion-report) is a tool for parsing outputs from fusion detection tools.
The score is explained here: <https://github.com/Clinical-Genomics/fusion-report/blob/master/docs/score.md>. Summary:

The weights for databases are as follows:

- COSMIC (50)
- MITELMAN (50)
- FusionGDB2 (0)

The final formula for calculating score is:

$$
score = 0.5 * \sum_{tool}^{tools} f(fusion, tool)*w(tool) + 0.5 * \sum_{db}^{dbs} g(fusion, db)*w(db)
$$

All tools have the same weight.

### Salmon

<details markdown="1">
<summary>Output files</summary>

- `salmon`
  - `<sample>`

</details>

Folder containing the quantification results

### Kallisto

<details markdown="1">
<summary>Output files</summary>

- `kallisto`
  - `<sample>.kallisto_quant.fusions.txt`

</details>

Quantifying abundances of transcripts from bulk and single-cell RNA-Seq data, or more generally of target sequences using high-throughput sequencing reads.

### Vcf_collect

<details markdown="1">
<summary>Output files</summary>

- `vcf_collect`
  - `<sample>_fusion_data.vcf` - contains the fusions in vcf format with collected statistics.

Vcf-collect takes as input the results of fusion-report and fusioninspector. That means fusions from all tools are aggregated. Fusioninspector applies a filter so it is possible some fusions detected by a caller are not filtered out by fusioninspector. In those cases, vcf-collect will display the fusions, but a lot of data will be missing as fusioninspector performs the analysis for each fusion.

</details>

[Megafusion](https://github.com/J35P312/MegaFusion) converts RNA fusion files to SV VCF and collects statistics and metrics in a VCF file.

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

### Picard

<details markdown="1">
<summary>Output files</summary>

Picard CollectRnaMetrics and picard MarkDuplicates share the same output directory.

- `picard`
  - `<sample>.MarkDuplicates.metrics.txt` - metrics from MarkDuplicates
  - `<sample>_rna_metrics.txt` - metrics from CollectRnaMetrics
  - `<sample>_insert_size_metrics.txt.txt` - metrics from CollectInsertSizeMetrics
  - `<sample>.bam` - BAM file with marked duplicates

</details>

### Samtools

#### Samtools sort

Samtools sort is used to sort BAM files from STAR_FOR_STARFUSION (for arriba visualisation)

<details markdown="1">
<summary>Output files</summary>

- `samtools_sort_for_<arriba>`
  - `<sample>(_chimeric)_sorted.bam` - sorted BAM file

</details>

#### Samtools index

Samtools index is used to index BAM files from STAR_FOR_ARRIBA (for arriba visualisation) and STAR_FOR_STARFUSION (for QC)

<details markdown="1">
<summary>Output files</summary>

- `samtools_for_<arriba/qc>`
  - `<sample>.(Aligned.sortedByCoord).out.bam.bai` -

</details>

### STAR

STAR is used to align to genome reference

STAR is run for 3 tools:

For `arriba` with the parameters:

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

For `STAR-fusion` with the parameters:

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
--chimOutType Junctions \
--quantMode GeneCounts
```

> STAR_FOR_STARFUSION uses `${params.ensembl}/Homo_sapiens.GRCh38.${params.ensembl_version}.chr.gtf` whereas STAR_FOR_ARRIBA uses `${params.ensembl_ref}/Homo_sapiens.GRCh38.${params.ensembl_version}.gtf`

<details markdown="1">
<summary>Output files</summary>

**Common**

- `star_for_<tool>`
- `<sample>.Log.final.out`
- `<sample>.Log.progress.out`
- `<sample>.SJ.out.tab`

**For arriba:**

- `<sample>.Aligned.out.bam`

  **For starfusion:**

- `<sample>.Aligned.sortedByCoord.out.bam`
- `<sample>.Chimeric.out.junction`
- `<sample>.ReadsPerGene.out.tab`

</details>

The STAR index is generated with `--sjdbOverhang ${params.read_length - 1}`, params.read_length default is 100.

### STAR-fusion

<details markdown="1">
<summary>Output files</summary>

- `starfusion`
  - `<sample>.starfusion.fusion_predictions.tsv` - contains the identified fusions
  - `<sample>.starfusion.abridged.tsv` - contains the identified fusions abridged
  - `starfusion.abridged.coding_effect.tsv`

</details>

### StringTie

<details markdown="1">
<summary>Output files</summary>

- `stringtie/<sample>/stringtie.merged.gtf` - merged gtf from annotation and stringtie output gtfs
</details>

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
