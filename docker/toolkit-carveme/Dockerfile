FROM quay.io/metagenomics/toolkit-python:3.13-0
RUN pip install carveme==1.6.5 pyscipopt \
	&& wget http://github.com/bbuchfink/diamond/releases/download/v2.1.13/diamond-linux64.tar.gz \
	&& tar xzvf diamond-linux64.tar.gz \
	&& chmod a+x diamond \
	&& mv diamond /usr/local/bin/

