FROM quay.io/biocontainers/virsorter:2.2.3--pyhdfd78af_1
LABEL maintainer="Jiarong Guo"
LABEL version="2.2.3"
RUN mamba install -y -c conda-forge wget
RUN virsorter setup -d /db && chmod -R 755 /db
VOLUME [ "/app" ]
WORKDIR /app
ENTRYPOINT [ "virsorter" ]
