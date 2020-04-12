FROM nfcore/base:1.9

LABEL authors="Martin Proks" \
    description="Docker image containing all requirements for nfcore/rnafusion pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/fusioncatcher_v1.20/bin:$PATH
