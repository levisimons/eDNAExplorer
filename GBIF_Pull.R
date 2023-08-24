#This script downloads the most recent copy of the GBIF database to a local drive.
require(gbifdb)
require(dplyr)
readRenviron(".env")
gbif_dir <- Sys.getenv("GBIF_HOME")
system(paste("sudo rm -rf ",gbif_dir,"/occurrence",sep=""))
gbif_download(dir=gbif_dir) #Run monthly to create GBIF mirror
#By default, this will download to the dir given by gbif_dir() and can be changed to a different path.
