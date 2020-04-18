FROM nfcore/base:1.9

LABEL authors="Martin Proks" \
    description="Docker image containing all requirements for nfcore/rnafusion pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/fusion-inspector_v1.3.1-star2.7.0f/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name fusion-inspector_v1.3.1-star2.7.0f > fusion-inspector_v1.3.1-star2.7.0f.yml