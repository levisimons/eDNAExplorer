r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos=r)
install.packages("crayon")
install.packages("BiocManager", repos = "http://cran.us.r-project.org")
install.packages("ade4")
BiocManager::install("Biostrings")
BiocManager::install("phyloseq")
