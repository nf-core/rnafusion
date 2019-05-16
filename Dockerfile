FROM nfcore/base
LABEL authors="Martin Proks <martin.proks@scilifelab.se>" \
      description="Docker image containing all requirements for nf-core/rnafusion pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-rnafusion-1.0.2/bin:$PATH
