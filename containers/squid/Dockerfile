FROM nfcore/base:1.9

LABEL authors="Martin Proks" \
    description="Docker image containing all requirements for nfcore/rnafusion pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-rnafusion-squid_1.5-star2.7.1a/bin:$PATH

RUN cd /opt/conda/envs/nf-core-rnafusion-squid_1.5-star2.7.1a/bin \
    && wget https://raw.githubusercontent.com/Kingsford-Group/squid/f45c9025d41cffd982ecbbdd52844e5a4f074de9/utils/AnnotateSQUIDOutput.py \
    && chmod +x /opt/conda/envs/nf-core-rnafusion-squid_1.5-star2.7.1a/bin/AnnotateSQUIDOutput.py \
    && ln -s /opt/conda/envs/nf-core-rnafusion-squid_1.5-star2.7.1a/bin/python3 /bin/python

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-rnafusion-squid_1.5-star2.7.1a > nf-core-rnafusion-squid_1.5-star2.7.1a.yml