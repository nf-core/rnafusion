FROM nfcore/base:1.12.1
LABEL authors="Martin Proks" \
      description="Docker image containing all software requirements for the nf-core/rnafusion pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-rnafusion-1.3.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-rnafusion-1.3.0dev > nf-core-rnafusion-1.3.0dev.yml

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron
