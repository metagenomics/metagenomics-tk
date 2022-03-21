FROM ubuntu:21.04
ENV VERSION=3.0.4
RUN ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime && \
    apt update && apt install wget git python3-pip -y && \
    git clone -b ${VERSION} https://github.com/phac-nml/mob-suite.git && \
    cd mob-suite && git checkout mob-3.0.3 && cd .. && \
    apt install libcurl4-openssl-dev libssl-dev -y && \
    pip3 install Cython numpy && \
    apt install mash ncbi-blast+ -y && \
    cd mob-suite && python3 setup.py install && cd .. && rm -rf mob-suite  && \
    apt clean
