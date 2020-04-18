FROM nfcore/base:1.7

LABEL authors="Martin Proks" \
    description="Docker image containing all requirements for nfcore/rnafusion pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/fusioncatcher_v1.20/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name fusioncatcher_v1.20 > fusioncatcher_v1.20.yml