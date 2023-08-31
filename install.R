if (!requireNamespace("phyloseq", quietly = TRUE)) {
    print("Attempting to install phyloseq")

    # Install dependencies
    deps <- c("Rhdf5lib", "rhdf5filters", "rhdf5", "biomformat")
    for (pkg in deps) {
        print(paste("Attempting to install", pkg))
        if (pkg == "Rhdf5lib") {
            install.packages("Rhdf5lib",
                configure.args = c("--with-hdf5=/usr/lib/x86_64-linux-gnu/hdf5/serial"),
                repos = BiocManager::repositories(),
                type = "source"
            )
        } else {
            BiocManager::install(pkg)
        }
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop(paste("Failed to install", pkg, ". Halting execution."))
        }
    }

    BiocManager::install("phyloseq")
    if (!requireNamespace("phyloseq", quietly = TRUE)) {
        stop("Failed to install phyloseq. Halting execution.")
    }
} else {
    print("phyloseq is installed")
}
install.packages("datadogr", repos = "https://cloud.r-project.org/")
if (!requireNamespace("datadogr", quietly = TRUE)) {
    stop("Failed to install datadogr. Halting execution.")
}
install.packages("ggVennDiagram", repos = "https://cloud.r-project.org/")
if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
    stop("Failed to install ggVennDiagram. Halting execution.")
}
