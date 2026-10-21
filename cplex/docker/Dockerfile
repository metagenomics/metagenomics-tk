# ---- stage 1: CPLEX installer ----
FROM python:3.7-slim-bookworm AS cplex

RUN apt-get update \
 && apt-get install -y --no-install-recommends unzip \
 && rm -rf /var/lib/apt/lists/*

COPY --from=eclipse-temurin:11-jre /opt/java/openjdk /opt/java/openjdk
ENV JAVA_HOME=/opt/java/openjdk
ENV PATH=/opt/java/openjdk/bin:${PATH}

COPY cos_installer-*.bin /tmp/installer
COPY install.properties /tmp/install.properties

RUN chmod u+x /tmp/installer \
 && /tmp/installer -i silent -f /tmp/install.properties \
 && test -d /opt/CPLEX/cplex/python/3.7 \
 || (echo "CPLEX install failed; contents of /opt:" && ls -la /opt && exit 1)


# ---- stage 2: libsbml 5.18.0 (the version sybilSBML 3.1.2 was written for) ----
FROM python:3.7-slim-bookworm AS libsbml

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      ca-certificates wget unzip build-essential cmake \
      libxml2-dev zlib1g-dev libbz2-dev \
 && rm -rf /var/lib/apt/lists/*

RUN set -e; \
    mkdir /tmp/src; \
    if wget --no-verbose -O /tmp/libsbml.tar.gz \
         "https://downloads.sourceforge.net/project/sbml/libsbml/5.18.0/experimental/source/libSBML-5.18.0-Source.tar.gz" \
       && tar -xzf /tmp/libsbml.tar.gz -C /tmp/src; then \
      echo "got source from SourceForge"; \
    else \
      echo "SourceForge failed, trying Zenodo"; \
      rm -rf /tmp/src/*; \
      wget --no-verbose -O /tmp/libsbml.zip \
        "https://zenodo.org/record/2645216/files/libSBML-5.18.0-Source.zip"; \
      unzip -q /tmp/libsbml.zip -d /tmp/src; \
    fi; \
    SRC="$(dirname "$(find /tmp/src -maxdepth 2 -name CMakeLists.txt | head -n1)")"; \
    echo "source dir: $SRC"; \
    cmake -S "$SRC" -B /tmp/build \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=/opt/libsbml \
      -DCMAKE_INSTALL_LIBDIR=lib \
      -DENABLE_FBC=ON \
      -DWITH_LIBXML=ON -DWITH_EXPAT=OFF \
      -DWITH_ZLIB=ON -DWITH_BZIP2=ON; \
    cmake --build /tmp/build -j"$(nproc)"; \
    cmake --install /tmp/build; \
    ls /opt/libsbml/lib /opt/libsbml/include/sbml/packages/fbc | head -20


# ---- stage 3: final image ----
FROM python:3.7-slim-bookworm

ARG COSDIR=/opt/CPLEX
ARG CPX_PYVERSION=3.7
ENV DEBIAN_FRONTEND=noninteractive
ENV CPX_PYVERSION=${CPX_PYVERSION}

ENV GAPSEQ_VERSION=v1.2
ENV CARVEME_VERSION=1.5.1
ENV SMETANA_VERSION=1.1.0
ENV MEMOTE_VERSION=0.14.0
ENV DIAMOND_VERSION=v2.0.13

WORKDIR /

COPY --from=cplex /opt/CPLEX /opt/CPLEX
COPY --from=libsbml /opt/libsbml /opt/libsbml
RUN echo /opt/libsbml/lib > /etc/ld.so.conf.d/libsbml.conf && ldconfig

# libsbml5-dev is intentionally not installed, so Debian's 5.19 headers
# can't shadow the 5.18 build in /opt/libsbml
RUN apt-get update && apt-get install -y --no-install-recommends \
      ca-certificates wget git bc parallel \
      ncbi-blast+ exonerate bedtools barrnap \
      r-base-core r-base-dev r-recommended \
      libglpk-dev libcurl4-openssl-dev libxml2-dev \
      libxml2 zlib1g libbz2-1.0 \
 && apt-mark manual libxml2 zlib1g libbz2-1.0

RUN R -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); install.packages(c("data.table", "stringr", "glpkAPI", "getopt", "doParallel", "foreach", "R.utils", "stringi", "CHNOSZ", "jsonlite"))' \
 && R -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Biostrings", ask = FALSE, update = FALSE)'

# glpkAPI was archived on CRAN (2026-03-20), so install it from the archive
RUN wget -q https://cran.r-project.org/src/contrib/Archive/glpkAPI/glpkAPI_1.3.4.1.tar.gz \
 && R CMD INSTALL glpkAPI_1.3.4.1.tar.gz \
 && rm -f glpkAPI_1.3.4.1.tar.gz \
 && R -q -e 'stopifnot(requireNamespace("glpkAPI"))'

RUN wget -q https://cran.r-project.org/src/contrib/Archive/sybil/sybil_2.2.0.tar.gz \
 && wget -q https://cran.r-project.org/src/contrib/Archive/sybilSBML/sybilSBML_3.1.2.tar.gz \
 && wget -q https://cran.r-project.org/src/contrib/Archive/cplexAPI/cplexAPI_1.4.0.tar.gz \
 && R CMD INSTALL sybil_2.2.0.tar.gz

RUN R CMD INSTALL --configure-args="--with-sbml-include=/opt/libsbml/include --with-sbml-lib=/opt/libsbml/lib" sybilSBML_3.1.2.tar.gz

RUN cd ${COSDIR}/cplex/python/${CPX_PYVERSION}/x86-64_linux \
 && python${CPX_PYVERSION} setup.py install

RUN R CMD INSTALL --configure-args=" \
      PKG_CFLAGS='-m64 -fPIC' \
      PKG_CPPFLAGS='-I${COSDIR}/cplex/include' \
      PKG_LIBS='-L${COSDIR}/cplex/lib/x86-64_linux/static_pic -lcplex -lm -pthread -ldl'" \
      cplexAPI_1.4.0.tar.gz

RUN R -e 'pk <- c("data.table","stringr","glpkAPI","getopt","doParallel","foreach","R.utils","stringi","CHNOSZ","jsonlite","Biostrings","sybil","sybilSBML","cplexAPI"); miss <- pk[!sapply(pk, requireNamespace, quietly=TRUE)]; if (length(miss)) stop("Missing R packages: ", paste(miss, collapse=", "))'

RUN rm -f sybil_2.2.0.tar.gz sybilSBML_3.1.2.tar.gz cplexAPI_1.4.0.tar.gz \
 && rm -rf \
      ${COSDIR}/concert ${COSDIR}/cpoptimizer ${COSDIR}/doc ${COSDIR}/opl \
      ${COSDIR}/python ${COSDIR}/Uninstall \
      ${COSDIR}/cplex/bin ${COSDIR}/cplex/examples ${COSDIR}/cplex/include \
      ${COSDIR}/cplex/lib ${COSDIR}/cplex/matlab ${COSDIR}/cplex/readmeUNIX.html \
 && (ls -d ${COSDIR}/cplex/python/* | grep -v "/${CPX_PYVERSION}$" | xargs -r rm -rf)

# ---- locales, gapseq, diamond, python tools ----
# apt-get update is needed here because the base image ships without package lists
RUN apt-get update \
 && apt-get install -y --no-install-recommends locales \
 && echo "LC_ALL=en_US.UTF-8" >> /etc/environment \
 && echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
 && echo "LANG=en_US.UTF-8" > /etc/locale.conf \
 && locale-gen en_US.UTF-8 \
 && dpkg-reconfigure locales \
 && git clone --depth 1 --branch ${GAPSEQ_VERSION} https://github.com/jotech/gapseq \
 && wget -q http://github.com/bbuchfink/diamond/releases/download/${DIAMOND_VERSION}/diamond-linux64.tar.gz \
 && tar xzf diamond-linux64.tar.gz && rm diamond-linux64.tar.gz \
 && pip install --no-cache-dir python-libsbml==5.20.2 swig \
 && pip install --no-cache-dir smetana==${SMETANA_VERSION} carveme==${CARVEME_VERSION} cobra==0.27 depinfo==1.7.0 memote==${MEMOTE_VERSION}

ENV PATH="${PATH}:/:/gapseq"

# final cleanup, after everything that needs the compiler and wget
RUN apt-get purge -y --auto-remove r-base-dev wget \
 && rm -rf /var/lib/apt/lists/*
