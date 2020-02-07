FROM nfcore/base

LABEL authors="rickard.hammaren@scilifelab.se, phil.ewels@scilifelab.se, martin.proks@scilifelab.se" \
    description="Docker image containing all requirements for nfcore/rnafusion pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

ENV PATH /opt/conda/envs/star-fusion_v1.6.0/bin:$PATH

# ctat-genome-lib-builder
RUN cd /opt/conda/envs/star-fusion_v1.6.0/ && \
    wget https://github.com/NCIP/ctat-genome-lib-builder/archive/2830cd708c5bb9353878ca98069427e83acdd625.zip && \
    unzip 2830cd708c5bb9353878ca98069427e83acdd625.zip && \
    cd ctat-genome-lib-builder-2830cd708c5bb9353878ca98069427e83acdd625 && \
    rm /opt/conda/envs/star-fusion_v1.6.0/bin/prep_genome_lib.pl

ENV PATH /opt/conda/envs/star-fusion_v1.6.0/ctat-genome-lib-builder-2830cd708c5bb9353878ca98069427e83acdd625/:$PATH

RUN ln -s /opt/conda/envs/star-fusion_v1.6.0/lib/libssl.so /opt/conda/envs/star-fusion_v1.6.0/lib/libssl.so.1.0.0 && \
    ln -s /opt/conda/envs/star-fusion_v1.6.0/lib/libcrypto.so /opt/conda/envs/star-fusion_v1.6.0/lib/libcrypto.so.1.0.0
