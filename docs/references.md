# nfcore/rnafusion: Download references for tools

Downloading references manually is a tedious long process. To make the pipeline easier to work with, we provide a script to download all necessary references for fusion detection tools.

> **TL;DR:** Make sure to download the correct references for your need!

```bash
nextflow run nf-core/rnafusion/download-references.nf --help
```

## Download all references

```bash
# Replace <COSMIC_USER> and <COSMIC_PASSWD> with yout credentials from COSMIC
nextflow run nf-core/rnafusion/download-references.nf \
  --download_all \
  --outdir <PATH> \
  --cosmic_usr <COSMIC_USER> --cosmic_passwd <COSMIC_PASSWD>
```

## Download specific references

```bash
# Example of downloading specific tools
nextflow run nf-core/rnafusion/download-references.nf \
--arriba \
--outdir <PATH>
```

## Tool reference requirements

| Tool             |        FASTA       |         GTF        |     STAR-index     |        Other       |
| ---------------- | :----------------: | :----------------: | :----------------: | :----------------: |
| Arriba           | :white_check_mark: | :white_check_mark: | :white_check_mark: | `custom_reference` |
| EricScript       |         :x:        |         :x:        |         :x:        | `custom_reference` |
| FusionCatcher    |         :x:        |         :x:        |         :x:        | `custom_reference` |
| Fusion-Inspector | :white_check_mark: | :white_check_mark: | :white_check_mark: |  `ctat_genome_lib` |
| fusion-report    |         :x:        |         :x:        |         :x:        |     `databases`    |
| Pizzly           |         :x:        | :white_check_mark: | :white_check_mark: |       `cDNA`       |
| Squid            |         :x:        | :white_check_mark: | :white_check_mark: |          -         |
| Star-Fusion      | :white_check_mark: | :white_check_mark: | :white_check_mark: |  `ctat_genome_lib` |
