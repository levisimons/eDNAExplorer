# eDNAExplorer
Script prototypes for the eDNA Explorer project

[![Run on Google Cloud](https://storage.googleapis.com/cloudrun/button.svg)](https://console.cloud.google.com/cloudshell/editor?shellonly=true&cloudshell_image=gcr.io/cloudrun/button&cloudshell_git_repo=https://github.com/maxogden/eDNAExplorer.git&revision=max/dockerize)

# Workflow
## Metadata
Users create a project and upload sample metadata.  This include sample locations, dates, and associated spatial uncertainties, along with any additional environmental variables captured in the field.

The full GBIF database is automatically and locally mirrored on the first day of each month using [this script](https://github.com/levisimons/eDNAExplorer/blob/main/GBIF_Pull.R)

The sample metadata is used to run a tool to extract additional environmental variables associated with each sampling location and date.  The metadata extractor is in development [here](https://github.com/MetadataExtractor).  This tool generates an extracted metadata file.

The sample and extracted metadata are merged using [this script](https://github.com/levisimons/eDNAExplorer/blob/main/eDNAExplorer_Metabarcoding_Metadata_Initializer.R).  This script is also used to determine which countries and states / provinces each sample lie within.  This information is used to enable comparisons between traditional observations of biodiversity and eDNA data.
