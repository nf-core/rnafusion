FROM nfcore/base:1.9

LABEL authors="Martin Proks" \
    description="Docker image containing all requirements for nfcore/rnafusion pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/pizzly_v0.37.3/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name pizzly_v0.37.3 > pizzly_v0.37.3.yml