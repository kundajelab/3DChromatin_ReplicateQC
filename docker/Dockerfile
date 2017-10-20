FROM ubuntu:16.04
MAINTAINER Kundaje Lab <oursu@stanford.edu>

RUN apt-get update --fix-missing
RUN apt-get install git -y
RUN apt-get install bzip2

RUN mkdir -p /src
WORKDIR /src

RUN git clone https://github.com/kundajelab/3DChromatin_ReplicateQC
WORKDIR /src/3DChromatin_ReplicateQC

ENV CONDA_DIR /opt/conda
ENV PATH $CONDA_DIR/bin:$PATH
RUN apt-get update --fix-missing
RUN apt-get install -y wget git libhdf5-dev g++ graphviz
RUN wget --quiet https://repo.continuum.io/miniconda/Miniconda2-4.1.11-Linux-x86_64.sh -O ~/miniconda.sh && /bin/bash ~/miniconda.sh -b -p $CONDA_DIR && rm ~/miniconda.sh
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh	
ENV PATH /opt/conda/bin:$PATH

ENV repodir /src/3DChromatin_ReplicateQC
RUN mkdir -p /src/3DChromatin_ReplicateQC/software
RUN git clone https://github.com/kundajelab/genomedisco $repodir/software/genomedisco

