# Use an official Ubuntu as a parent image
FROM 576661951794.dkr.ecr.us-west-1.amazonaws.com/reports:latest

# Copy R dependencies file
COPY install.dev.R /tmp/install.dev.R
RUN conda run -n reports /bin/bash -c "Rscript /tmp/install.dev.R"