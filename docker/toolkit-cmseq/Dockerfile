FROM ubuntu:18.04
WORKDIR /usr/local/bin
ENV LC_ALL=C.UTF-8 
RUN apt update &&  \
     apt install -y wget zlib1g-dev libbz2-dev liblzma-dev python3 python3-pip && \ 
     wget https://github.com/SegataLab/cmseq/archive/1.0.tar.gz && \ 
     tar xzvf 1.0.tar.gz  --strip 1 && \
     pip3 install Cython && pip3 install  . -vv
