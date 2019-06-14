FROM continuumio/miniconda3
ENV VERSION 0.0.1

RUN useradd --create-home sim
WORKDIR /home/sim

RUN conda install -c bioconda numpy=1.16.4 viennarna=2.4.13 && conda clean -a
RUN git clone https://github.com/desiro/silentMutations.git 

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y inkscape && \
    apt-get install -y default-jre && \
    apt-get clean

WORKDIR /home/sim/silentMutations
ENV PATH=${PATH}:"/home/sim/silentMutations/"
ENV PATH=${PATH}:"/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/"

USER sim
