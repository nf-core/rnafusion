# ![rnafusion](https://raw.githubusercontent.com/nf-core/rnafusion/master/docs/images/rnafusion_logo.png.png)

## *NB* This is work in progress - use at your own risk!

Attempt at setting up a pipeline for fusion detection as an addon to the previously existing nf-core/rnaseq pipeline.

## Usage:
### ` --reads`
Path to input reads
**NB: Must be enclosed in quotes!**
```
 --reads 'path/to/data/sample_*_{1,2}.fastq'
```

### `--star-fusion`
If this flag is set to true then the pipeline will run STAR-Fusion.
True by default

### `--fusioncatcher`

If this flag is set to true then the pipeline will run Fusion Catcher.
True by default

More than one of these can be true at a time.

Example usage if you wish to run both STAR-fusion and Fusion Catcher   
```
nextflow run nf-core/rnafusion/main.nf --star --fusioncatcher --reads 'SRR5665649_{1,2}*.gz'
```

### `--fc_extra_options`

Invoking this flag allows for extra options to be sent directly to Fusion Catcher.
