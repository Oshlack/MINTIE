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

library(dplyr)
library(data.table)

options(stringsAsFactors = FALSE,
        error = function(){dump.frames("compare_eq_classes_debug", to.file = TRUE); q()})

# load helper methods
args <- commandArgs(trailingOnly=FALSE)
file.arg <- grep("--file=", args, value=TRUE)
incl.path <- gsub("--file=(.*)compare_eq_classes.R", "\\1de_methods.R", file.arg)
source(incl.path, chdir=TRUE)

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 3) {
    print("Invalid arguments.\n\n")
    args <- c("--help")
}

if("--help" %in% args) {
    cat("
        Load in Salmon equivalence classes from two samples,
        match them and perform differential expression analysis
        between equivalence class counts between case and control
        samples.

        Usage:
        Rscript compare_eq_classes.R <case_name> <ec_matrix> <tx_ref_fasta> <output> --FDR=<value> --minCPM=<value> --minLogFC=<value> --test)\n\n

        All flags are optonal. The defaults are:
        FDR = 0.05
        minCPM = 0.1
        minLogFC = 5\n\n
    ")

    q(save="no")
}

MIN_CONTROLS <- 2
REC_CONTROLS <- 10

GENE_REGEX <- "gene_symbol:([a-zA-Z0-9.-]+)"

FDR <- 0.05
minCPM <- 0.1
minLogFC <- 2

#############################################################
# Parse arguments
#############################################################

args <- commandArgs(trailingOnly=TRUE)
case_name <- args[1]
ec_matrix_file <- args[2]
tx_ref_fasta <- args[3]
outfile <- args[4]
test_mode <- length(grep("--test", args)) > 0

# optional flags
set_arg <- function(argname) {
    flag <- paste0("--", argname, "=")
    arg <- grep(flag, args, value=T)
    if(!length(arg) == 0) {
        value <- as.numeric(strsplit(arg, "=")[[1]][2])
        if(is.na(value)) {
            stop(paste("Invalid", argname, "value."))
        }
        assign(argname, value, envir = .GlobalEnv)
    }
}
set_arg("FDR")
set_arg("minCPM")
set_arg("minLogFC")

#############################################################
# load data
#############################################################

print("Reading input data...")
ec_matrix <- fread(ec_matrix_file, sep="\t")
n_controls <- ncol(ec_matrix) - 4
if(!test_mode) {
    if(n_controls < MIN_CONTROLS) {
        stop(paste("Insufficient controls. Please run MINTIE with at least", MIN_CONTROLS, "controls."))
    } else if(n_controls < REC_CONTROLS) {
        print(paste("WARNING: you are running MINTIE with fewer than", REC_CONTROLS, "controls.",
                    "Adding more control samples will improve results."))
    }
}

# get reference transcript IDs
txs <- fread(tx_ref_fasta, header=FALSE, sep="\n")
txs <- txs[grep("^>", txs$V1),]
txs <- sapply(txs$V1, function(x){strsplit(x, " ")[[1]][1]})
txs <- as.character(sapply(txs, gsub, pattern=">", replacement=""))

#############################################################
# Prepare data for DE analysis
#############################################################

print("Extracting ECs associated with only novel contigs...")
tx_ec <- data.table(distinct(ec_matrix[,c("ec_names", "transcript")]))
tx_ec$novel <- !tx_ec$transcript%in%txs
novel_contig_ecs <- tx_ec[, all(novel), by="ec_names"]
novel_contig_ecs <- unique(novel_contig_ecs$ec_names[novel_contig_ecs$V1])
if(length(novel_contig_ecs) == 0) {
    stop("No novel ECs were found. Pipeline cannot continue.")
}

# make collapsed EC > contig lookup table
tx_ec <- tx_ec[tx_ec$ec_names%in%novel_contig_ecs,]
tx_ec <- tx_ec[, paste(transcript, collapse=":"), by="ec_names"]
colnames(tx_ec)[2] <- "contigs"

#############################################################
# Perform DE
#############################################################

print("Performing differential expression analysis...")
de_results <- run_edgeR(case_name, ec_matrix, tx_ec, dirname(outfile),
                        cpm_cutoff=minCPM, qval=FDR, min_logfc=minLogFC, test=test_mode)
if(nrow(de_results) == 0) {
    stop("Invalid output result obtained from differential expression. Please check all input data.")
}

#############################################################
# Compile and write results
#############################################################

print("Compiling and writing results...")

# separate out contigs
contigs <- sapply(de_results$contigs, strsplit, split=":")
de_results <- de_results[rep(1:nrow(de_results), as.numeric(sapply(contigs, length))),]
de_results$contig <- as.character(unlist(contigs))

write.table(de_results, outfile, row.names=F, quote=F, sep="\t")
