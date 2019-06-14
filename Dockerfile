FROM continuumio/miniconda3
ENV VERSION 0.0.1

RUN useradd --create-home sim
WORKDIR /home/sim

RUN conda install -c bioconda numpy=1.16.4 viennarna=2.4.13 && conda clean -a
RUN git clone https://github.com/desiro/silentMutations.git 

RUN apt-get update && apt-get upgrade
RUN apt-get install inkscape=0.92.4
RUN apt-get install oracle-java11-installer

WORKDIR /home/sim/silentMutations
ENV PATH=${PATH}:"/home/sim/silentMutations"

USER sim
