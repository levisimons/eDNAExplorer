# Use an official Ubuntu as a parent image
FROM ubuntu:latest

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
  dirmngr

# Add CRAN Repository for R 4.3.x
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" -y
RUN apt-get install -y r-base r-base-dev

# Install AWS CLI
RUN wget "https://d1vvhvl2y92vvt.cloudfront.net/awscli-exe-linux-x86_64.zip" -O awscliv2.zip && \
  unzip awscliv2.zip && \
  ./aws/install -i /root/aws
ENV PATH="/root/aws/bin:${PATH}"

# # Copy env.yaml and install.R to the image
COPY install.R /tmp/install.R
COPY install_biocmanager.R /tmp/install_biocmanager.R

# Install R packages such as rhdf5 and Phyloseq
RUN /bin/bash -c "Rscript /tmp/install.R"
RUN /bin/bash -c "Rscript /tmp/install_biocmanager.R"

# Set the working directory
WORKDIR /project

# Copy entrypoint script into the image
COPY ./scripts/entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

# Set the entry point
ENTRYPOINT ["/entrypoint.sh"]



