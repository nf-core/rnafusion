# nfcore/rnafusion: Tools

## Using tool
```bash
cd utils && sh download.sh
```

## Manual download

### STAR-Fusion

- [Homepage](https://github.com/STAR-Fusion/STAR-Fusion/wiki)
- Minimum requirement: ~30GB RAM
```bash
wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz -O GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz
tar -xvzf GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz
```

Update the config file to include the directory
```groovy
params {
  star_fusion_ref = "/path/to/GRCh38_v27_CTAT_lib_Feb092018/ctat_genome_lib_build_dir"
}
```

## Fusioncatcher

- [Homepage](https://github.com/ndaniel/fusioncatcher/blob/master/doc/manual.md)
- Minimum requirement: ~60GB RAM
```bash
wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v90.tar.gz.aa
wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v90.tar.gz.ab
wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v90.tar.gz.ac
wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v90.tar.gz.ad
cat human_v90.tar.gz.* | tar xz
```

Update the config file to include the directory
```groovy
params {
  fusioncatcher_ref = '/path/to/human_v90'
}
```
## Ericscript

- [Homepage](https://sites.google.com/site/bioericscript/getting-started)
- Genome: [https://drive.google.com/file/d/0B9s__vuJPvIiUGt1SnFMZFg4TlE/view]
```bash
wget https://raw.githubusercontent.com/circulosmeos/gdown.pl/master/gdown.pl \
&& chmod +x gdown.pl \
&& ./gdown.pl "https://drive.google.com/uc?export=download&confirm=qgOc&id=0B9s__vuJPvIiUGt1SnFMZFg4TlE" ericscript_db_homosapiens_ensembl84.tar.bz2 \
&& rm gdown.pl
```
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
wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz \

# annotation
wget ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz && gunzip Homo_sapiens.GRCh38.94.gtf.gz
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
```bash
mkdir -p igenomes/Homo_sapiens/NCBI/GRCh38/
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/NCBI/GRCh38/ .
```
## FusionInspector

- [Homepage](https://github.com/FusionInspector/FusionInspector/wiki)
- Uses reference genome from STAR-Fusion (ctat_genome_lib_build_dir)

## Custom summary report

The final summary report is made so that the user can customize it easily. These are the supported changes:
- Custom title, date format
- Add institution logo
- Change styling of the report

```js
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