# eDNAExplorer

Script prototypes for the eDNA Explorer project

[![Run on Google Cloud](https://storage.googleapis.com/cloudrun/button.svg)](https://console.cloud.google.com/cloudshell/editor?shellonly=true&cloudshell_image=gcr.io/cloudrun/button&cloudshell_git_repo=https://github.com/maxogden/eDNAExplorer.git&revision=max/dockerize)

# Development

## Configure Local Environment Variables

Run the following to setup your local environment file:

`cp .env.sample .env`

Edit the `.env` file with the actual credentials necessary to connect to AWS, postgres, and redis.

## Developing Locally with a Pre-Built Container

You can develop locally now with our pre-built container:
https://hub.docker.com/repository/docker/jimjeffers/edna-explorer/general

For your convenience we've supplied some shell scripts to easily execute commands or bash into the container to test things. First and foremost ensure you've configured your local environment variables. For the time being we will be using production environment variables until a staging environment is setup. So please use the `.env` from the jet stream instance so that you can connect to the jetstream bucket and the supabase database.

### Running Local Scripts

You can quickly and easily run scripts via the `./scripts/run.sh` shell script. Here's an example:

```
./scripts/run.sh Rscript --vanilla eDNAExplorer_Alpha_Metabarcoding.R "clk06lckc0001jn0f9hw7rb86" "2020-07-24" "2021-02-20" "CO1_Metazoa" "25" "species" "1000" "0.00003" "None" "ghm" "Chao1"
```

### Run an Interactive Shell

If you need to test something within the environment to ensure things are working as expected or to experiment installing different packages without rebuilding the docker container you can use the `shell.sh` script as follows:

```
./scripts/shell.sh`
```

## Build the Docker Image Locally

Provided you have [docker installed](https://docs.docker.com/desktop/) - you can build a local instance of the conda environment with phyloseq pre-installed with the following command:

`docker build . --platform linux/amd64 -t edna-explorer`

After successfully building the image you can run an interactive shell to execute local scripts via:

`docker run --platform linux/amd64 -v $PWD:/project -it edna-explorer /bin/bash`

From there you can simply `cd ./project` and execute R scripts directly via the command line:

```
cd ./project
Rscript --vanilla eDNAExplorer_Alpha_Metabarcoding.R "clk3duk760001ml0fug8cvmp0" "2020-07-24" "2021-02-20" "CO1_Metazoa" "25" "species" "1000" "0.00003" "None" "ghm" "Chao1"
```

# Workflow

## Metadata

Users create a project and upload sample metadata. This include sample locations, dates, and associated spatial uncertainties, along with any additional environmental variables captured in the field.

The full GBIF database is automatically and locally mirrored on the first day of each month using [this script](https://github.com/levisimons/eDNAExplorer/blob/main/GBIF_Pull.R)

The sample metadata is used to run a tool to extract additional environmental variables associated with each sampling location and date. The metadata extractor is in development [here](https://github.com/MetadataExtractor). This tool generates an extracted metadata file.

### Metabarcoded samples

The sample and extracted metadata are merged and stored in a PostgreSQL database using [this script](https://github.com/levisimons/eDNAExplorer/blob/main/eDNAExplorer_Metabarcoding_Metadata_Initializer.R). This script is also used to determine which countries and states / provinces each sample lie within. This information is used to enable comparisons between traditional observations of biodiversity and eDNA data.

### qPCR samples

(Lorem impsum and then some)

## Taxonomic data
