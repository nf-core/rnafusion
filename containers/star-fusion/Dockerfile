FROM nfcore/base:1.9

LABEL authors="Martin Proks" \
    description="Docker image containing all requirements for nfcore/rnafusion pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-rnafusion-star-fusion_1.8.1/bin:$PATH

# ctat-genome-lib-builder
RUN cd /opt/conda/envs/nf-core-rnafusion-star-fusion_1.8.1/ && \
    wget https://github.com/NCIP/ctat-genome-lib-builder/archive/2830cd708c5bb9353878ca98069427e83acdd625.zip && \
    unzip 2830cd708c5bb9353878ca98069427e83acdd625.zip && \
    cd ctat-genome-lib-builder-2830cd708c5bb9353878ca98069427e83acdd625 && \
    rm /opt/conda/envs/nf-core-rnafusion-star-fusion_1.8.1/bin/prep_genome_lib.pl

ENV PATH /opt/conda/envs/nf-core-rnafusion-star-fusion_1.8.1/ctat-genome-lib-builder-2830cd708c5bb9353878ca98069427e83acdd625/:$PATH

RUN ln -s /opt/conda/envs/nf-core-rnafusion-star-fusion_1.8.1/lib/libssl.so /opt/conda/envs/nf-core-rnafusion-star-fusion_1.8.1/lib/libssl.so.1.0.0 && \
    ln -s /opt/conda/envs/nf-core-rnafusion-star-fusion_1.8.1/lib/libcrypto.so /opt/conda/envs/nf-core-rnafusion-star-fusion_1.8.1/lib/libcrypto.so.1.0.0

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-rnafusion-star-fusion_1.8.1 > nf-core-rnafusion-star-fusion_1.8.1.yml