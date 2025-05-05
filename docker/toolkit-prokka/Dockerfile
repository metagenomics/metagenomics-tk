FROM ubuntu:18.04

# metadata
LABEL website="https://github.com/tseemann/prokka"
LABEL license="https://github.com/tseemann/prokka#licence"

# install dependencies
RUN apt-get update && \
    apt-get -y --no-install-recommends install libdatetime-perl libxml-simple-perl libdigest-md5-perl git default-jre bioperl hmmer wget \
                       zlib1g-dev python liblzma-dev libbz2-dev xz-utils curl g++ make \
                       less bzip2 gzip libidn11-dev &&\ 
    apt-get clean && apt-get autoclean && rm -rf /var/lib/apt/lists/*

#get Bioperl
#RUN cpan Bio::Perl

# install bedtools 2.29.0 since >=2.27.0 is required for barrnap and the apt-get package is 2.25.0
# dependencies required for bedtools: zlib1g-dev python liblzma-dev libbz2-dev xz-utils curl g++
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.29.0/bedtools-2.29.0.tar.gz && \
  tar -zxf bedtools-2.29.0.tar.gz && \
  rm bedtools-2.29.0.tar.gz && \
  cd bedtools2 && \
  make

# add bedtools to PATH for barrnap test
ENV PATH="$PATH:/bedtools2/bin"

# install barrnap
RUN wget https://github.com/tseemann/barrnap/archive/0.9.tar.gz &&\
    tar -zxf 0.9.tar.gz &&\
    rm 0.9.tar.gz &&\
    cd barrnap-0.9 &&\
    make test

#get prokka
RUN git clone https://github.com/metagenomics/prokka

# install tbl2asn manually since the one included with prokka is expired.
RUN mv /prokka/tbl2asn /usr/bin/tbl2asn && \
    chmod 755 /usr/bin/tbl2asn

# put prokka in PATH & 
# set PATH to prioritize /usr/bin/tbl2asn before /prokka/binaries/tbl2asn
ENV PATH="$PATH:/prokka/bin"

#setupdb
RUN prokka --setupdb

WORKDIR /data

