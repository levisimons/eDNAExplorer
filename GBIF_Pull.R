# This script downloads the most recent copy of the GBIF database to a local drive.
require(gbifdb)
require(dplyr)
library(sentryR)
readRenviron(".env")

configure_sentry(
    dsn = Sys.getenv("SENTRY_DSN"),
    app_name = "r-report-service", app_version = "1.1.0",
    environment = Sys.getenv("APP_ENV"),
    runtime = NULL
)
cat("Configured Sentry: ", Sys.getenv("SENTRY_DSN"), "\n")

gbif_dir <- Sys.getenv("GBIF_HOME")
cat("GBIF Home: ", gbif_dir, "\n")
cat("Removing old GBIF data\n")
system(paste("rm -rf ", gbif_dir, "/occurrence", sep = ""))
cat("Downloading new GBIF data\n")
gbif_download(dir = gbif_dir) # Run monthly to create GBIF mirror
