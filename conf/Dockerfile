FROM openjdk:8

LABEL authors="rickard.hammaren@scilifelab.se, phil.ewels@scilifelab.se" \
    description="Docker image containing all requirements for NGI-RNAfusion pipeline"



#Fusioncatcher 
RUN wget http://sf.net/projects/fusioncatcher/files/bootstrap.py -O bootstrap.py && python bootstrap.py -t --download -y


#FusionInspector
RUN wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.4.0.tar.gz \
    make \
    make plugins



