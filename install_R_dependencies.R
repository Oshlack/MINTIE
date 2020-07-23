repos="http://cran.r-project.org"
if (!require("dplyr")) {
    install.packages("dplyr", repos=repos)
}
if (!require("data.table")) {
    install.packages("data.table", repos=repos)
}
if (!require("readr")) {
    install.packages("readr", repos=repos)
}
if (!require("jsonlite")) {
    install.packages("jsonlite", repos=repos)
}
if (!require("tximport") | !require("edgeR")) {
    r_version = paste(R.Version()$major, strsplit(R.Version()$minor, '\\.')[[1]][1], sep='.')
    if(as.numeric(r_version) < 3.5) {
        source("https://bioconductor.org/biocLite.R")
        biocLite("tximport")
        biocLite("edgeR")
    } else {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager")
        }
        BiocManager::install("tximport")
        BiocManager::install("edgeR")
    }
}
