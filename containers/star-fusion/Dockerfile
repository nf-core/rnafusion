FROM nfcore/base:1.9

LABEL authors="Martin Proks" \
    description="Docker image containing all requirements for nfcore/rnafusion pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-rnafusion-star-fusion_1.8.1/bin:$PATH

# FusionInspector
ENV PATH /opt/conda/envs/nf-core-rnafusion-star-fusion_1.8.1/lib/STAR-Fusion/FusionInspector:$PATH

# ctat-genome-lib-builder
ENV PATH /opt/conda/envs/nf-core-rnafusion-star-fusion_1.8.1/lib/STAR-Fusion/ctat-genome-lib-builder:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-rnafusion-star-fusion_1.8.1 > nf-core-rnafusion-star-fusion_1.8.1.yml