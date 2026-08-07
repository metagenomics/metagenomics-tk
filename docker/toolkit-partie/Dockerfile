FROM ubuntu:22.10

# Install tzdata to set the timezone to UTC (required for R)
# A non-interactive frontend is used to avoid tzdata asking for user input
RUN apt-get update &&\
    DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata

RUN apt-get update && \
    apt-get install -y r-base perl curl git wget && \
    R -e "install.packages('randomForest', repos='http://cran.rstudio.com/')" && \
    wget https://github.com/gmarcais/Jellyfish/releases/download/v2.3.0/jellyfish-2.3.0.tar.gz && \
    tar -xvzf jellyfish-2.3.0.tar.gz && \
    cd jellyfish-2.3.0 && \
    ./configure && make && make install && \
    cd .. && rm -rf jellyfish-2.3.0* && \
    apt-get install -y libjellyfish-2.0-2 && \
    wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.4/bowtie2-2.4.4-linux-x86_64.zip && \
    unzip bowtie2-2.4.4-linux-x86_64.zip && \
    mv bowtie2-2.4.4-linux-x86_64/bowtie2* /usr/local/bin/ && \
    git clone https://github.com/linsalrob/partie.git && \
    wget https://github.com/lh3/seqtk/archive/refs/tags/v1.3.tar.gz && \
    tar -xvzf v1.3.tar.gz && \
    cd seqtk-1.3 && make && cp seqtk /usr/local/bin/ && cd .. && rm -rf seqtk-1.3* && \
    wget http://moo.nac.uci.edu/~hjm/zotkill.pl -O /usr/local/bin/zotkill.pl && chmod +x /usr/local/bin/zotkill.pl

COPY partie.pl /partie/partie.pl

ENV PATH="/usr/local/bin:/partie:${PATH}"
#
WORKDIR /partie