#This script downloads the most recent copy of the GBIF database to a local drive.
require(gbifdb)
require(dplyr)
gbif_download() #Run monthly to create GBIF mirror
#By default, this will download to the dir given by gbif_dir() and can be changed to a different path.
