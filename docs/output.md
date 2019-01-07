# nf-core/rnafusion: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.


## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [FastQC](#fastqc) - read quality control
* [Fusioncatcher](#fusioncatcher)
* [Ericscript](#ericscript)
* [MultiQC](#multiqc) - aggregate report, describing results of the whole pipeline

## FastQC
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the `trim_galore` directory.

**Output directory: `results/fastqc`**

* `sample_fastqc.html`
  * FastQC report, containing quality metrics for your untrimmed raw fastq files
* `zips/sample_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images

## Star-Fusion

**Output directory: `results/tools/StarFusion`**
* `star-fusion.fusion_predictions.tsv`
  * contains all the predicted gene fusions

For more info check the [documentation](https://github.com/STAR-Fusion/STAR-Fusion/wiki#Outputs).

## Fusioncatcher

**Output directory: `results/tools/Fusioncatcher`**
* `final-list_candidate-fusion-genes.txt`
  * contains all the predicted gene fusions

For more info check the [documentation](https://github.com/ndaniel/fusioncatcher/blob/master/doc/manual.md#62---output-data-output-data).

## Ericscript

**Output directory: `results/tools/Ericscript/tmp`**
* `fusions.results.total.tsv`
  * contains all the predicted gene fusions

For more info check the [documentation](https://sites.google.com/site/bioericscript/getting-started).

## Pizzly

**Output directory: `results/tools/Pizzly`**
* `pizzly_fusions.unfiltered.json`
  * contains all the predicted gene fusions

For more info check the [documentation](https://github.com/pmelsted/pizzly#output).

## Squid

**Output directory: `results/tools/Squid`**
* `fusions_annotated.txt`
  * contains all the predicted gene fusions

For more info check the [documentation](https://github.com/Kingsford-Group/squid#output-specification).

## Fusion Inspector

**Output directory: `results/tools/FusionInspector`**
* `finspector.fa`
  * the candidate fusion-gene contigs (if you copy things elsewhere, make sure to als copy the index file finspector.fa.fai)
* `finspector.bed`
  * the reference gene structure annotations for fusion partners
* `finspector.junction_reads.bam`
  * alignments of the breakpoint-junction supporting reads.
* `finspector.spanning_reads.bam`
  * alignments of the breakpoint-spanning paired-end reads.

To visualize fusion genes in [IGV tool](https://software.broadinstitute.org/software/igv/igvtools) first create a genome `Menu->Genomes->Create .genome File`, choose name and description, then choose the following files:

* `finspector.fa`
  * Make sure the index file finspector.fa.fai is in the same folder
* `finspector.gtf`
  * use this for 'Genes'
* `cytoBand.txt`
  * use this for 'optional Cytoband'

Add the bam files by choosing `File->Load from File` and make sure to select your generated mini genome in the upper-left corner. 

For more info and help check [wiki page](https://github.com/FusionInspector/FusionInspector/wiki).

## MultiQC
[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `results/multiqc`**

* `Project_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `Project_multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see http://multiqc.info
