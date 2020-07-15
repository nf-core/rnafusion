FROM nfcore/base:1.7

LABEL authors="Martin Proks" \
    description="Docker image containing all requirements for nfcore/rnafusion pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-rnafusion-fusioncatcher_1.20/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-rnafusion-fusioncatcher_1.20 > nf-core-rnafusion-fusioncatcher_1.20.yml