FROM python:3.7

ARG COSDIR=/opt/CPLEX

ARG CPX_PYVERSION=3.7

COPY cos_installer-*.bin /tmp/installer
COPY install.properties /tmp/install.properties
RUN  chmod u+x /tmp/installer && apt-get update && apt-get install -y default-jre \
  && /tmp/installer -f /tmp/install.properties \
  && apt install -y ncbi-blast+ git r-base-core exonerate libglpk-dev bedtools barrnap bc parallel libcurl4-openssl-dev libsbml5-dev  \
  && R -e 'install.packages(c("data.table", "stringr", "glpkAPI", "sybil", "getopt", "doParallel", "foreach", "R.utils", "stringi", "CHNOSZ", "jsonlite"))' \
  && R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Biostrings")' \
  && wget https://cran.r-project.org/src/contrib/Archive/sybilSBML/sybilSBML_3.1.2.tar.gz &&  R CMD INSTALL sybilSBML_3.1.2.tar.gz \ 
  && wget https://cran.r-project.org/src/contrib/Archive/cplexAPI/cplexAPI_1.4.0.tar.gz \
  && apt-get remove -y --purge default-jre && apt-get -y --purge autoremove \
  && cd ${COSDIR}/cplex/python/${CPX_PYVERSION}/x86-64_linux \
  &&      python${CPX_PYVERSION} setup.py install && cd - \ 
  && R CMD INSTALL --configure-args=" \
    PKG_CFLAGS='-m64 -fPIC' \
    PKG_CPPFLAGS='-I/opt/CPLEX/cplex/include' \
    PKG_LIBS='-L/opt/CPLEX/cplex/lib/x86-64_linux/static_pic -lcplex -lm -pthread -ldl'" cplexAPI_1.4.0.tar.gz \
  && rm -f /tmp/installer /tmp/install.properties \
  && rm -rf \
   ${COSDIR}/concert \
   ${COSDIR}/cpoptimizer \
   ${COSDIR}/doc \
   ${COSDIR}/opl \
   ${COSDIR}/python \
   ${COSDIR}/Uninstall \
   ${COSDIR}/cplex/bin \
   ${COSDIR}/cplex/examples \
   ${COSDIR}/cplex/include \
   ${COSDIR}/cplex/lib \
   ${COSDIR}/cplex/matlab \
   ${COSDIR}/cplex/readmeUNIX.html \
  && ls -d ${COSDIR}/cplex/python/* | grep -v ${CPX_PYVERSION} | xargs rm -rf


ENV CPX_PYVERSION ${CPX_PYVERSION}

RUN DEBIAN_FRONTEND=noninteractive apt-get install -y locales \
    && echo "LC_ALL=en_US.UTF-8" >> /etc/environment \
    && echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
    && echo "LANG=en_US.UTF-8" > /etc/locale.conf \
    && locale-gen en_US.UTF-8 \
  && DEBIAN_FRONTEND=noninteractive dpkg-reconfigure locales && git clone https://github.com/jotech/gapseq \
  && wget http://github.com/bbuchfink/diamond/releases/download/v2.0.13/diamond-linux64.tar.gz  \
  && tar xzf diamond-linux64.tar.gz \
  && pip install smetana carveme memote


ENV PATH="${PATH}:/gapseq"
