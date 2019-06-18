FROM continuumio/miniconda3
ENV VERSION 0.0.1

RUN useradd --create-home sim
WORKDIR /home/sim

RUN conda install -c bioconda numpy=1.16.4 viennarna=2.4.13 && conda clean -a
RUN git clone https://github.com/desiro/silentMutations.git 

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y openjdk-8-jdk && \
    apt-get install -y ant && \
    apt-get clean

# Fix certificate issues
RUN apt-get update && \
    apt-get install ca-certificates-java && \
    apt-get clean && \
    update-ca-certificates -f;

# Setup JAVA_HOME -- useful for docker commandline
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN export JAVA_HOME

WORKDIR /home/sim/silentMutations
ENV PATH=${PATH}:"/home/sim/silentMutations"
ENV PATH=${PATH}:"/home/sim/silentMutations/VARNAv3-93.jar"

USER sim
