install_and_verify_package <- function(
    package_name,
    ncpus = 8) {
    r <- getOption("repos")
    r["CRAN"] <- "http://cran.us.r-project.org"
    r["cloud"] <- "https://cloud.r-project.org/"
    r["duck"] <- "https://duckdb.r-universe.dev"
    options(repos = r)
    install.packages(package_name, Ncpus = ncpus)
    if (!requireNamespace(package_name, quietly = TRUE)) {
        stop(paste("Failed to install", package_name, ". Halting execution."))
    }
}

# Usage
install_and_verify_package("datadogr")
install_and_verify_package("ggVennDiagram")
install_and_verify_package("plotly")
install_and_verify_package("tidyr")
install_and_verify_package("dplyr")
install_and_verify_package("lubridate")
install_and_verify_package("jsonlite")
install_and_verify_package("RPostgreSQL")
install_and_verify_package("anytime")
install_and_verify_package("devtools")
install_and_verify_package("duckdb")
install_and_verify_package("duckdbfs")
install_and_verify_package("sentryR")
devtools::install_github("ropensci/gbifdb")
devtools::install_github("cboettig/minioclient")
