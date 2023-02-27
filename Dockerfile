FROM rocker/r-base

RUN apt-get update -qq && apt-get install -y \
  git-core \
  libssl-dev \
  libxml2-dev \
  libcurl4-gnutls-dev \ 
  libsodium-dev

RUN install2.r plumber
COPY ["./install.R", "./install.R"]
RUN ["Rscript", "./install.R"]
COPY [".", "./"]

ENTRYPOINT ["R", "-e", "pr <- plumber::plumb(commandArgs()[4]); pr$run(host='0.0.0.0', port=as.numeric(Sys.getenv('PORT')))"]
CMD ["api.R"]
