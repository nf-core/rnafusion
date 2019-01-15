#!/bin/bash
if [ $# -eq 0 ]; then
    echo "No path specified"
    exit 1
else
    DEST=$1
    if [ -d "$DEST" ]; then
        mkdir $DEST/rnafusion_references
        cd $DEST/rnafusion_references

        # STAR-Fusion
        echo 'Downloading STAR-Fusion references'
        mkdir star_fusion_ref && cd star_fusion_ref \
        && wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz -O GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz \
        && tar -xvzf GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz && rm GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz

        # Fusioncatcher
        echo 'Downloading Fusioncatcher references'
        cd $cwd \
        && mkdir fusioncatcher_ref && cd fusioncatcher_ref \
        && wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v90.tar.gz.aa \
        && wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v90.tar.gz.ab \
        && wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v90.tar.gz.ac \
        && wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v90.tar.gz.ad \
        && cat human_v90.tar.gz.* | tar xz \
        && rm human_v90.tar*

        # Ericscript
        echo 'Downloading Ericscript references'
        cd $cwd \
        && mkdir ericscript_ref && cd ericscript_ref \
        && wget https://raw.githubusercontent.com/circulosmeos/gdown.pl/dfd6dc910a38a42d550397bb5c2335be2c4bcf54/gdown.pl \
        && chmod +x gdown.pl \
        && ./gdown.pl "https://drive.google.com/uc?export=download&confirm=qgOc&id=0B9s__vuJPvIiUGt1SnFMZFg4TlE" ericscript_db_homosapiens_ensembl84.tar.bz2 \
        && tar jxf ericscript_db_homosapiens_ensembl84.tar.bz2 \
        && rm gdown.pl ericscript_db_homosapiens_ensembl84.tar.bz2

        # Pizzly
        echo 'Downloading Pizzly references'
        cd $cwd \
        && mkdir pizzly_ref && cd pizzly_ref \
        && wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz \
        && wget ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz && gunzip Homo_sapiens.GRCh38.94.gtf.gz

        # iGenome references
        # mkdir -p $cwd/igenomes/Homo_sapiens/NCBI/GRCh38/
        # aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/NCBI/GRCh38/ .
    else
        echo "Specified path doesn't exists"
    fi
fi