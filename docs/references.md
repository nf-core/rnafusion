# nfcore/rnafusion: Download references for tools

## Use nextflow helper script

Downloading references manually is a tedious long process. To make the pipeline easier to work with, we provide a script to download all necessary references for fusion detection tools.

> **TL;DR:** Make sure to download the correct references for your need!

```bash
nextflow run nf-core/rnafusion/download-references.nf --help
```

```bash
# Example to download all tools
nextflow run nf-core/rnafusion/download-references.nf \
  -profile <PROFILE> \
  --arriba \
  --ericscript \
  --fusioncatcher \
  --fusion_report --cosmic_usr <COSMIC_USER> --cosmic_passwd <COSMIC_PASSWD> \
  --pizzly \
  --star_fusion \
  --outdir <PATH>
```

```bash
# Example of downloading specific tools
nextflow run nf-core/rnafusion/download-references.nf \
  -profile <PROFILE> \
  --arriba \
  --fusion_report --cosmic_usr <COSMIC_USER> --cosmic_passwd <COSMIC_PASSWD> \
  --outdir <PATH>
```

## Tool reference requirements

### STAR-Fusion

- STAR-index
- FASTA
- GTF
- ctat_genome_lib_build_dir (built using `prep_genome_lib.pl`)

### Fusioncatcher

- Custom reference (human)
- FASTA
- GTF

### Ericscript

- Custom reference
- FASTA
- GTF

### Pizzly

- Transcript
- FASTA
- GTF

### Squid

- STAR-index
- GTF file

### Fusion-Inspector

- ctat_genome_lib_build_dir
- FASTA
- GTF

### fusion-report

- internal databases
