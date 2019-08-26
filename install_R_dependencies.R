repos="http://cran.r-project.org"
if (!require("dplyr")) {
    install.packages("dplyr", repos=repos)
}
if (!require("seqinr")) {
    install.packages("seqinr", repos=repos)
}
if (!require("data.table")) {
    install.packages("data.table", repos=repos)
}
if (!require("edgeR")) {
    r_version = paste(R.Version()$major, strsplit(R.Version()$minor, '\\.')[[1]][1], sep='.')
    if(as.numeric(r_version) < 3.5) {
        source("https://bioconductor.org/biocLite.R")
        biocLite("edgeR")
    } else {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager")
        }
        BiocManager::install("edgeR")
    }
}
