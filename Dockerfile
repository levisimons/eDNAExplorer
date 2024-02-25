# Use an official Ubuntu as a parent image
FROM ubuntu:jammy

# Install tzdata as it usually blocks the installation of other packages due to prompt
RUN apt-get clean && \
  apt-get update --fix-missing && \
  apt-get upgrade -y
RUN DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata

# Run package updates and install packages
RUN apt-get install -y \
  wget \
  bzip2 \
  unzip \
  libfontconfig1-dev \
  libgdal-dev \
  libcurl4-gnutls-dev \
  libssl-dev \
  libxml2-dev \
  libudunits2-dev \
  zlib1g-dev \
  libproj-dev \
  libcairo2-dev \
  libv8-dev \
  libsodium-dev \
  libpq-dev \
  libsqlite3-dev \
  libboost-dev \
  libgeos-dev \
  software-properties-common \
  dirmngr \
  awscli \
  parallel \
  git

# Add CRAN Repository for R 4.3.x
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" -y
RUN apt-get install -y r-base r-base-dev

# Install conda
RUN wget -P /tmp/ "https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh" && \
  bash "/tmp/Miniconda3-py38_4.12.0-Linux-x86_64.sh" -b -p /usr/local/miniconda

COPY env.yml /tmp/env.yml
ENV PATH="/usr/local/miniconda/bin:$PATH"

# Copy git files
COPY . /home/ubuntu/eDNAExplorer/

# Create reports env
RUN conda env create -f /tmp/env.yml -n reports && \
  conda init

RUN echo "source activate reports" >> /root/.bashrc

# Copy R dependencies file
COPY install.R /tmp/install.R
COPY install_biocmanager.R /tmp/install_biocmanager.R

# Install R packages such as rhdf5 and Phyloseq
RUN conda run -n reports /bin/bash -c "Rscript /tmp/install.R"
RUN conda run -n reports /bin/bash -c "Rscript /tmp/install_biocmanager.R"

# Install minio client for GBIF Pull
RUN wget https://dl.min.io/client/mc/release/linux-amd64/mc -O /usr/local/bin/mc && \
  chmod +x /usr/local/bin/mc

# Create a symlink for mc binary to be accessible as expected by minioclient
RUN mkdir -p /root/.local/share/R/minioclient/ && \
  ln -s /usr/local/bin/mc /root/.local/share/R/minioclient/mc

# Set the working directory
WORKDIR /home/ubuntu/eDNAExplorer