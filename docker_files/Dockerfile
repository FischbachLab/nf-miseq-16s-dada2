# Dockerfile - https://github.com/blekhmanlab/docker_dada2/blob/master/Dockerfile

# Base Image
FROM bioconductor/bioconductor_docker:RELEASE_3_12

LABEL authors="Xiandong Meng" \
      description="16s purity check using DADA2 " \
      software.name="DADA2" \
      software.description="16s" \
      software.website="https://benjjneb.github.io/dada2/" \
      container.category="tool"

USER root:root

# Install base utilities
#RUN apt-get update && \
#    apt-get install -y build-essentials  && \
#    apt-get install -y wget &&
#    apt-get clean && \
#    rm -rf /var/lib/apt/lists/*

# install iniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
     /bin/bash ~/miniconda.sh -b -p /opt/conda

ENV PATH "/opt/conda/bin:$PATH"

# this "version" label reflects the installed version of DADA2:
LABEL org.label-schema.version = "1.18.0"

RUN R -e 'BiocManager::install("dada2", version="3.12", update=TRUE, ask=FALSE)'
RUN R -e 'BiocManager::install("phyloseq", version="3.12", update=TRUE, ask=FALSE)'
#RUN R -e 'BiocManager::install("decontam", version="3.12", update=TRUE, ask=FALSE)'


# Update conda to latest version.
RUN conda update -n base -c defaults conda

# Add conda channels
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# Install software packages
#RUN conda install -c bioconda -y bioconductor-dada2 \
  #&& conda install -c bioconda -y bioconductor-phyloseq  \
RUN conda install -y numpy pandas python \
   && conda install -c conda-forge -y awscli \
   && conda clean -afy

# Install time and memory usage tool (timem)
RUN conda install -c jrmadsen -y timemory

ADD fix_summary.py /usr/local/bin/fix_summary.py
RUN chmod +x /usr/local/bin/fix_summary.py
ADD 16s_wrapper.sh /usr/local/bin/16s_wrapper.sh
RUN chmod +x /usr/local/bin/16s_wrapper.sh
ADD 16S.R /usr/local/bin/16S.R
RUN chmod +x /usr/local/bin/16S.R

RUN mkdir -p /work
WORKDIR /work

# Get Repo
COPY . .

RUN chmod -R +rx ./

CMD ["pwd"]
