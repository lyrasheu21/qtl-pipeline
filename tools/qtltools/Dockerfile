# QTLtools v1.3.1 built from source, for running the pipeline's QTLtools rules
# in a container:
#   docker build -t qtltools:1.3.1 -f tools/qtltools/Dockerfile .
#   # config:  qtl: {qtltools_binary: QTLtools, container: "docker://<registry>/qtltools:1.3.1"}
#   snakemake --profile profiles/local --sdm conda apptainer ...
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential gfortran curl ca-certificates bzip2 gzip \
      libboost-iostreams-dev libboost-program-options-dev libgsl-dev \
      zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev \
      libpcre2-dev \
    && rm -rf /var/lib/apt/lists/*

COPY tools/qtltools/install_qtltools.sh /opt/install_qtltools.sh
RUN bash /opt/install_qtltools.sh /opt/qtltools \
    && cp /opt/qtltools/bin/QTLtools /usr/local/bin/QTLtools \
    && rm -rf /opt/qtltools/build

ENTRYPOINT []
CMD ["QTLtools", "--help"]
