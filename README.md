# ![NGI-RNAfusion](https://raw.githubusercontent.com/SciLifeLab/NGI-RNAfusion/master/docs/images/NGI-RNAfusion_logo.png)

## *NB* This is work in progress - use at your own risk!

Attempt at setting up a pipeline for fusion detection as an addon to the previously existing NGI-RNAseq pipeline.

## Usage:
### ` --reads`
Path to input reads 
**NB: Must be enclosed in quotes!**
```
 --reads 'path/to/data/sample_*_{1,2}.fastq'
```

### `--star`
If this flag is set to true then the pipeline will run STAR-Fusion.

### `--fusioncatcher`

If this flag is set to true then the pipeline will run Fusion Catcher.

More than one of these can be true at a time. 

Example usage if you wish to run both STAR-fusion and Fusion Catcher   
```
nextflow run /NGI-RNAfusion/main.nf --star --fusioncatcher --reads 'SRR5665649_{1,2}*.gz'
```

### `--fc_extra_options`

Invoking this flag allows for extra options to be sent directly to Fusion Catcher. This is implemented since Fusion Catcher has too many potential options too be having their own flag in this pipeline. 

---

[![SciLifeLab](https://raw.githubusercontent.com/SciLifeLab/NGI-RNAfusion/master/docs/images/SciLifeLab_logo.png)](http://www.scilifelab.se/)
[![National Genomics Infrastructure](https://raw.githubusercontent.com/SciLifeLab/NGI-RNAfusion/master/docs/images/NGI_logo.png)](https://ngisweden.scilifelab.se/)

---
