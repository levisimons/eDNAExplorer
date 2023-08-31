# Use an official Ubuntu as a parent image
FROM ubuntu:latest

# Run package updates and install packages
RUN apt-get clean && \
  apt-get update --fix-missing && \
  apt-get upgrade -y && \
  apt-get install -y \
  wget \
  bzip2 \
  unzip \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev \
  libudunits2-dev \
  zlib1g-dev \
  libhdf5-dev

# Install AWS CLI
RUN wget "https://d1vvhvl2y92vvt.cloudfront.net/awscli-exe-linux-x86_64.zip" -O awscliv2.zip && \
  unzip awscliv2.zip && \
  ./aws/install -i /root/aws
ENV PATH="/root/aws/bin:${PATH}"

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
  bash ~/miniconda.sh -b -p /opt/conda && \
  rm ~/miniconda.sh

# Add conda to PATH
ENV PATH /opt/conda/bin:$PATH

# Copy env.yaml and install.R to the image
COPY env.yml /tmp/env.yml
COPY install.R /tmp/install.R

# Create a Conda environment from env.yml
RUN conda env create -f /tmp/env.yml

# Activate the Conda environment and run install.R
RUN echo "source activate $(head -1 /tmp/env.yml | cut -d' ' -f2)" > ~/.bashrc
ENV PATH /opt/conda/envs/$(head -1 /tmp/env.yml | cut -d' ' -f2)/bin:$PATH

# Install R packages such as rhdf5 and Phyloseq
RUN /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && conda activate GBIF_env && Rscript /tmp/install.R"

# Set the working directory
WORKDIR /project

# Copy entrypoint script into the image
COPY ./scripts/entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

# Set the entry point
ENTRYPOINT ["/entrypoint.sh"]



