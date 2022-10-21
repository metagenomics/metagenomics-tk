FROM continuumio/miniconda3

# Create the environment:
RUN git clone https://github.com/metagentools/MetaCoAG.git && cd MetaCoAG && git checkout 34c16943f3cd0a9c942711af6992ede7bec36e1e

RUN conda env create -f MetaCoAG/environment.yml

# Make RUN commands use the new environment:
RUN echo "conda activate metacoag" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

# The code to run when container is started:
ENTRYPOINT ["conda", "activate", "metacoag"]

ENV PATH="/MetaCoAG:${PATH}"
