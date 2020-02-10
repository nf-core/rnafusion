FROM nfcore/base:1.8
LABEL authors="Martin Proks" \
      description="Docker image containing all software requirements for the nf-core/rnafusion pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-rnafusion-dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-rnafusion-dev > nf-core-rnafusion-dev.yml
