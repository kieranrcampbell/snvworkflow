FROM bioconductor/release_core2

RUN apt-get update && \
    apt-get install -y curl


# Install miniconda to /miniconda
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh
RUN bash Miniconda-latest-Linux-x86_64.sh -p /miniconda -b
RUN rm Miniconda-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda

# Create environment
RUN conda create -n py27 python=2.7 anaconda

RUN apt-get update && apt install -y python3-pip && \
    pip3 install snakemake	  

# Source environment
RUN /bin/bash -c "source /miniconda/bin/activate py27"

# Install necessary packages
RUN conda install -c anaconda pandas && \
    conda install -c bioconda pysam && \
    conda install -c anaconda scipy && \
    conda install -c anaconda numpy
    
# Install R packages required
RUN Rscript -e "install.packages('devtools', dependencies=TRUE)" && \
    Rscript -e "install.packages('tidyverse', dependencies=TRUE)" && \
    Rscript -e "devtools::install_github('jeremystan/aargh')"
  