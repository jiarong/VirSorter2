FROM quay.io/biocontainers/virsorter:latest
LABEL maintainer="Jiarong Guo"
LABEL version="2.2.3"
RUN virsorter config --init-source --db-dir /work/projects/Cyverse/iVirus/VS2/db
#RUN virsorter setup --db /db
VOLUME [ "/app" ]
WORKDIR /app
ENTRYPOINT [ " virsorter " ]
