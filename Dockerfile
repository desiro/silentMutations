FROM continuumio/miniconda3
ENV VERSION 0.0.1

RUN useradd --create-home sim
WORKDIR /home/sim

RUN conda install -c bioconda numpy=1.16.4 viennarna=2.4.13 && conda clean -a
RUN conda install -c lb_arrakistx varna=3.93 && conda clean -a
RUN git clone https://github.com/desiro/silentMutations.git 

WORKDIR /home/sim/silentMutations
ENV PATH=${PATH}:"/home/sim/silentMutations"

#USER sim
