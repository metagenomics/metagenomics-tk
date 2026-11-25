FROM ubuntu:18.10
WORKDIR /sans
RUN apt update && apt install -y wget g++ git make python3.8 python3-dendropy && \
    git clone https://gitlab.ub.uni-bielefeld.de/gi/sans.git . && \
    chmod -R a+rw . && \
    make
