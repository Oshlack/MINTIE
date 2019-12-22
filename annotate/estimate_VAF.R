args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0 | (length(args) < 4 & !"--help" %in% args)) {
    cat("
        Invalid arguments.
    ")
    args <- c("--help")
}

if("--help" %in% args) {
    cat("
        Calculate rough VAF estimation from salmon quantification results

        Usage:
        Rscript estimate_VAF.R <quant.sf> <contig_info> <tx2gene> <outfile>\n\n
    ")

    q(save="no")
}

if (!require("tximport")) {
    r_version = paste(R.Version()$major, strsplit(R.Version()$minor, '\\.')[[1]][1], sep='.')
    if(as.numeric(r_version) < 3.5) {
        source("https://bioconductor.org/biocLite.R")
        biocLite("tximport")
    } else {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager")
        }
        BiocManager::install("tximport")
    }
}

library(tximport)
library(data.table)
library(dplyr)

options(stringsAsFactors = FALSE,
        error = function(){dump.frames("estimate_vaf", to.file = TRUE); q()})

#############################################################
# Parse arguments
#############################################################

quant_file <- args[1]
contig_info <- args[2]
tx2gene <- args[3]
outfile <- args[4]

#############################################################
# Load data
#############################################################

print("Loading data...")
cinfo <- read.delim(file = contig_info)

txi <- tximport(quant_file, type='salmon', countsFromAbundance = 'lengthScaledTPM', txOut=TRUE)
quant <- data.frame(txi$abundance)
colnames(quant) <- 'TPM'; quant$tx <- rownames(quant)

tx2g <- read.delim(tx2gene, header = FALSE)
colnames(tx2g) <- c('tx', 'gene')

#############################################################
# Calculate VAFs
#############################################################

print("Estimating VAFs...")
# make lookup table for which genes map to contigs
c2g <- data.frame(contig_id=cinfo$contig_id, gene=cinfo$overlapping_genes)
split_genes <- sapply(c2g$gene, function(x){strsplit(x, '\\||:')})
c2g <- data.frame(contig_id=rep(c2g$contig_id, sapply(split_genes, length)),
                  gene=unlist(split_genes))
c2g <- distinct(c2g)
c2g <- c2g[c2g$gene!='',]

# calculate wild-type TPMs (sum all WT transcripts)
wt_tpm <- merge(quant, tx2g, by='tx')
wt_tpm <- data.table(wt_tpm)[, sum(TPM), by='gene']
colnames(wt_tpm)[2] <- 'WT_TPM'

# merge to get relevant fields
x <- quant[!quant$tx%in%tx2g$tx & quant$TPM>0,]
x$contig_id <- x$tx
x <- merge(x, c2g, by='contig_id')
x <- merge(x, wt_tpm, by='gene')
x <- distinct(x)
x$tx <- NULL

# for contigs overlapping multiple genes, take average TPM of all WT genes
mean_tpm <- data.table(x)[, mean(WT_TPM), by='contig_id']
colnames(mean_tpm)[2] <- 'mean_WT_TPM'
x <- merge(x, mean_tpm, by='contig_id')

if (nrow(x) > 0) {
    # calculate VAF and write output
    x$VAF <- x$TPM / (x$TPM + x$mean_WT_TPM)
    write.table(x, file=outfile, row.names=FALSE, quote=FALSE, sep='\t')
} else {
    print("ERROR: no variants to output. Please check your tx2gene.txt reference file.")
}
