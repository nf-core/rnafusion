FROM nfcore/base:1.9

LABEL authors="Martin Proks" \
    description="Docker image containing all requirements for nfcore/rnafusion pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/arriba_v1.1.0/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name arriba_v1.1.0 > arriba_v1.1.0.yml