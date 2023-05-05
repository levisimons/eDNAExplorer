# eDNAExplorer
Script prototypes for the eDNA Explorer project

[![Run on Google Cloud](https://storage.googleapis.com/cloudrun/button.svg)](https://console.cloud.google.com/cloudshell/editor?shellonly=true&cloudshell_image=gcr.io/cloudrun/button&cloudshell_git_repo=https://github.com/maxogden/eDNAExplorer.git&revision=max/dockerize)

# Workflow
Users create a project and upload sample metadata.  This include sample locations, dates, and associated spatial uncertainties, along with any additional environmental variables captured in the field.
The full GBIF database is locally mirror on the first day of each month using [this script](https://github.com/levisimons/eDNAExplorer/blob/main/GBIF_Pull.R)
