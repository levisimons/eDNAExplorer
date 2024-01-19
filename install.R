install_and_verify_package <- function(
    package_name,
    repos_url = "https://cloud.r-project.org/",
    ncpus = 8) {
    install.packages(package_name, repos = repos_url, Ncpus = ncpus)
    if (!requireNamespace(package_name, quietly = TRUE)) {
        stop(paste("Failed to install", package_name, ". Halting execution."))
    }
}

# Usage
install_and_verify_package("datadogr")
install_and_verify_package("ggVennDiagram")
install_and_verify_package("ggplot2")
install_and_verify_package("plotly")
install_and_verify_package("tidyr")
install_and_verify_package("dplyr")
install_and_verify_package("lubridate")
install_and_verify_package("jsonlite")
install_and_verify_package("data.table")
install_and_verify_package("uuid")
install_and_verify_package("RPostgreSQL")
install_and_verify_package("DBI")
install_and_verify_package("anytime")
install_and_verify_package("digest")
install.packages('duckdb', repos=c('https://duckdb.r-universe.dev', 'https://cloud.r-project.org'))
install.packages('duckdbfs', repos=c('https://duckdb.r-universe.dev', 'https://cloud.r-project.org'))
devtools::install_github("ropensci/gbifdb")
devtools::install_github("cboettig/minioclient")

