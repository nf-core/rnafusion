# nfcore/rnafusion: Tools

## STAR-Fusion

- [Homepage](https://github.com/STAR-Fusion/STAR-Fusion/wiki)
- Minimum requirement: ~30GB RAM
```bash
wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz
tar -xvzf GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz
```

Update the config file to include the directory
```groovy
params {
  star_fusion_ref = "/path/to/ctat_genome_lib_build_dir"
}
```

## Fusioncatcher

- [Homepage](https://github.com/ndaniel/fusioncatcher/blob/master/doc/manual.md)
- Minimum requirement: ~60GB RAM
```bash
wget http://sf.net/projects/fusioncatcher/files/bootstrap.py -O bootstrap.py && python bootstrap.py --download
```

Update the config file to include the directory
```groovy
params {
  fusioncatcher_dir = '/path/to/fusioncatcher-ref'
}
```
## Ericscript

- [Homepage](https://sites.google.com/site/bioericscript/getting-started)
- Genome: [https://drive.google.com/file/d/0B9s__vuJPvIiUGt1SnFMZFg4TlE/view](https://drive.google.com/file/d/0B9s__vuJPvIiUGt1SnFMZFg4TlE/view)

Update the config file to include the directory
```groovy
params {
  ericscript_ref = '/path/to/ericscript_db_homosapiens_ensembl84'
}
```

## Pizzly

- [Homepage](https://github.com/pmelsted/pizzly)
```bash
# transcriptome
wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -O pizzly_ref/
# annotation
wget ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz && tar -xvzf Homo_sapiens.GRCh38.94.gtf.gz
```

Update the config file to include the directory
```groovy
params {
  pizzly_ref = '/path/to/pizzly-ref'
}
```

## Squid

- [Homepage](https://github.com/Kingsford-Group/squid)
- Requires STAR alignment, GTF, FASTA file

## FusionInspector

- [Homepage](https://github.com/FusionInspector/FusionInspector/wiki)
- Uses reference genome from STAR-Fusion (ctat_genome_lib_build_dir)

## Custom summary report

The final summary report is made so that the user can customize it easily. These are the supported changes:
- Custom title, date format
- Add institution logo
- Change styling of the report

```json
# Report customization configuration
report_title: 'Some really cool title'
institution: 'assets/img/scilifelab.jpg'
date_format: '%Y-%m-%d'
assets:
  css:
    - assets/css/something.css
  js:
    - assets/js/something.js
```