FROM nfcore/base:1.9

LABEL authors="Martin Proks" \
    description="Docker image containing all requirements for nfcore/rnafusion pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/squid_v1.5-star2.7.0f/bin:$PATH

RUN cd /opt/conda/envs/squid_v1.5-star2.7.0f/bin \
    && wget https://raw.githubusercontent.com/Kingsford-Group/squid/f45c9025d41cffd982ecbbdd52844e5a4f074de9/utils/AnnotateSQUIDOutput.py \
    && chmod +x /opt/conda/envs/squid_v1.5-star2.7.0f/bin/AnnotateSQUIDOutput.py \
    && ln -s /opt/conda/envs/squid_v1.5-star2.7.0f/bin/python3 /bin/python

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name squid_v1.5-star2.7.0f > squid_v1.5-star2.7.0f.yml