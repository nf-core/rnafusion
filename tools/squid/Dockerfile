FROM nfcore/base

LABEL authors="rickard.hammaren@scilifelab.se, phil.ewels@scilifelab.se, martin.proks@scilifelab.se" \
    description="Docker image containing all requirements for nfcore/rnafusion pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/squid_v1.5-star2.7.0f/bin:$PATH

RUN wget https://raw.githubusercontent.com/Kingsford-Group/squid/master/utils/AnnotateSQUIDOutput.py -O /AnnotateSQUIDOutput.py \
    && chmod +x /AnnotateSQUIDOutput.py \
    && ln -s /AnnotateSQUIDOutput.py /usr/local/bin \
    && ln -s /opt/conda/envs/squid_v1.5-star2.7.0f/bin/python3 /bin/python 

RUN ln -s /lib/x86_64-linux-gnu/libcrypt.so.1 /lib/x86_64-linux-gnu/libcrypto.so.1.0.0
