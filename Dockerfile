FROM ubuntu:17.10

WORKDIR /

RUN apt-get update && \
	apt-get upgrade -y && \
	apt-get install -y wget vim bzip2 less

RUN wget --quiet https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
	bash ~/miniconda.sh -b -p /opt/conda && \
	rm ~/miniconda.sh

ENV PATH /opt/conda/bin:$PATH

RUN conda update conda -y




WORKDIR /

COPY SpliceJunctionDiscovery.py /
COPY SpliceJunctionNormalization.py /
COPY SpliceJunctionFiltering.py /
COPY reference /reference/

RUN mkdir outputs

RUN apt-get update && apt-get install -y \
  bzip2 \
  g++ \
  libbz2-dev \
  liblzma-dev \
  make \
  ncurses-dev \
  wget \
  zlib1g-dev


ENV SAMTOOLS_INSTALL_DIR=/opt/samtools

WORKDIR /tmp
RUN wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2 && \
  tar --bzip2 -xf samtools-1.7.tar.bz2


WORKDIR /tmp/samtools-1.7
RUN ./configure --enable-plugins --prefix=$SAMTOOLS_INSTALL_DIR && \
  make all all-htslib && \
  make install install-htslib

WORKDIR /
RUN ln -s $SAMTOOLS_INSTALL_DIR/bin/samtools /usr/bin/samtools && \
  rm -rf /tmp/samtools-1.7

WORKDIR /

RUN pip install pandas