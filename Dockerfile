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
  libcurl4-openssl-dev \
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
  postgresql \ 
  postgresql-contrib \
  software-properties-common \
  dirmngr \
  awscli \
  parallel


# Add CRAN Repository for R 4.3.x
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" -y
RUN apt-get install -y r-base r-base-dev

# Install conda
RUN wget -P /tmp/ "https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh" && \
    bash "/tmp/Miniconda3-py38_4.12.0-Linux-x86_64.sh" -b -p /usr/local/miniconda

COPY env.yml /tmp/env.yml
COPY env_metadata.yml /tmp/env_metadata.yml
ENV PATH="/usr/local/miniconda/bin:$PATH"

# Create GBIF_env and metadata
RUN conda env create -f /tmp/env.yml -n GBIF_env && \
    conda env create -f /tmp/env_metadata.yml -n metadata && \
    conda init && \
    . /root/.bashrc

# Copy env.yaml and install.R to the image
COPY install.R /tmp/install.R
COPY install_biocmanager.R /tmp/install_biocmanager.R

# Install R packages such as rhdf5 and Phyloseq
RUN conda run -n GBIF_env /bin/bash -c "Rscript /tmp/install.R"
RUN conda run -n GBIF_env /bin/bash -c "Rscript /tmp/install_biocmanager.R"

# Set the working directory
WORKDIR /project

# Make port 8000 available to the world outside this container
EXPOSE 8000

# Activate the environment, and run the application
CMD ["bash", "-c", "source activate metadata && exec gunicorn --workers 3 --bind 0.0.0.0:8000 wsgi:app"]
