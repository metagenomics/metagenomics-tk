FROM quay.io/biocontainers/bwa-mem2:2.2.1--hd03093a_2 AS bwa-mem2

FROM ubuntu:20.04
ENV SAMTOOLS=samtools-1.12
RUN apt update && apt install -y bc libgomp1 libgetopt-complete-perl wget autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev
RUN wget -qO- https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2 | tar jxvf -  &&\
    cd ${SAMTOOLS}  &&\
    ./configure &&\
    make &&\
    make install &&\
    cd .. && rm -rf ${SAMTOOLS}
COPY --from=bwa-mem2 /usr/local/bin/bwa-mem2*    /usr/local/bin/ 
