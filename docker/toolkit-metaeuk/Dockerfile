FROM ubuntu:25.04
RUN apt update && apt install -y zstd wget
ENV RELEASE=7-bba0d80

RUN wget https://github.com/soedinglab/metaeuk/releases/download/${RELEASE}/metaeuk-linux-sse2.tar.gz \
	&& tar xzvf metaeuk-linux-sse2.tar.gz \
        && mv metaeuk/bin/metaeuk /usr/local/bin/metaeuk_sse2

RUN wget https://github.com/soedinglab/metaeuk/releases/download/${RELEASE}/metaeuk-linux-avx2.tar.gz \
	&& tar xzvf metaeuk-linux-avx2.tar.gz \
        && mv metaeuk/bin/metaeuk /usr/local/bin/metaeuk_avx2

RUN wget https://github.com/soedinglab/metaeuk/releases/download/${RELEASE}/metaeuk-linux-sse41.tar.gz  \
	&& tar xzvf metaeuk-linux-sse41.tar.gz \
        && mv metaeuk/bin/metaeuk /usr/local/bin/metaeuk_sse41
