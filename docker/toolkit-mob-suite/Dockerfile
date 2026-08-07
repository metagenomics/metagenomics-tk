FROM ubuntu:22.10
ENV VERSION=3.1.0
RUN ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime && \
    apt update && apt install software-properties-common  wget git -y && \
    add-apt-repository ppa:deadsnakes/ppa -y && \
    apt update && \
    apt install python3.8 python3-pip python3.8-distutils python3.8-dev -y  && \
    git clone -b ${VERSION} https://github.com/phac-nml/mob-suite.git && \
    apt install libcurl4-openssl-dev libssl-dev -y && \
    pip3 install Cython numpy && \
    apt install mash ncbi-blast+ -y && \
    cd mob-suite && python3.8 setup.py install && cd .. && rm -rf mob-suite  && \
    apt clean
