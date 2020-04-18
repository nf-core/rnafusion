FROM nfcore/base:1.9

LABEL authors="Martin Proks" \
    description="Docker image containing all requirements for nfcore/rnafusion pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/ericscript_v0.5.5/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name ericscript_v0.5.5 > ericscript_v0.5.5.yml