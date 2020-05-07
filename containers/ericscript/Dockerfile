FROM nfcore/base:1.9

LABEL authors="Martin Proks" \
    description="Docker image containing all requirements for nfcore/rnafusion pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-rnafusion-ericscript_0.5.5/bin:$PATH

# Ignore database check (https://github.com/nf-core/rnafusion/issues/119)
RUN echo 1 > /opt/conda/envs/nf-core-rnafusion-ericscript_0.5.5/share/ericscript-0.5.5-4/lib/data/_resources/.flag.dbexists

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-rnafusion-ericscript_0.5.5 > nf-core-rnafusion-ericscript_0.5.5.yml