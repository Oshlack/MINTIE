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
        Rscript estimate_VAF.R <ec_matrix_file> <quant_file> <contig_info_file> <tx_ref_fasta> <tx2gene_file> <outfile>\n\n
    ")

    q(save="no")
}

if (!require("tximport")) {
    r_version = paste(R.Version()$major, strsplit(R.Version()$minor, "\\.")[[1]][1], sep=".")
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

ec_matrix_file <- args[1]
quant_file <- args[2]
contig_info_file <- args[3]
tx_ref_fasta <- args[4]
tx2gene_file <- args[5]
outfile <- args[6]

#############################################################
# Load data
#############################################################

print("Loading data...")
txi <- tximport(quant_file, type="salmon", countsFromAbundance = "lengthScaledTPM", txOut=TRUE)
ec_matrix <- fread(ec_matrix_file)
cinfo <- fread(contig_info_file)
tx2g <- fread(tx2gene_file, col.names=c("transcript", "gene"))
txs <- fread(tx_ref_fasta, header=FALSE)

#############################################################
# Prepare data
#############################################################

print("Preparing data...")
# get overlapping gene info for novel annotated contigs
c2g <- data.frame(contig_id=cinfo$contig_id, gene=cinfo$overlapping_genes)
split_genes <- sapply(c2g$gene, function(x){strsplit(x, "\\||:")})
c2g <- data.frame(transcript=rep(c2g$contig_id, sapply(split_genes, length)),
                  gene=unlist(split_genes))
c2g <- distinct(c2g)
c2g <- c2g[c2g$gene!="",]

# match non-novel contig ECs to genes and combine with novel contigs
# we do this so that we have a gene mapping for every contig and transcript
like_refseq <- ec_matrix$transcript%like%"hg38_ncbiRefSeq"
if(any(like_refseq)) {
    ec_matrix[like_refseq, 'transcript'] <- sapply(ec_matrix$transcript[like_refseq],
                                                   function(x){strsplit(gsub("hg38_ncbiRefSeq_", "", x), "\\.")[[1]][1]})
}
ec2tx <- distinct(ec_matrix[,c("transcript", "ec_names")])
ec2g <- inner_join(ec2tx, tx2g, by="transcript")
if(nrow(ec2g) == 0) {
    stop("ERROR: no transcripts in the tx2gene reference match the supplied EC matrix! Please double check your reference and matrix file.")
}
tx2g <- distinct(ec2g[,c("transcript", "gene")])
tx2g <- rbind(tx2g, c2g)

# get list of all reference transcripts
# we have to get these from the fasta as some wildtype transcripts
# may not be in the tx2gene reference because they lack gene names
txs <- txs[grep("^>", txs$V1),]
txs <- sapply(txs$V1, function(x){strsplit(x, " ")[[1]][1]})
txs <- as.character(sapply(txs, gsub, pattern=">", replacement=""))
like_refseq <- txs%like%"hg38_ncbiRefSeq"
if(any(like_refseq)) {
    txs <- as.character(sapply(txs, function(x){strsplit(gsub("hg38_ncbiRefSeq_", "", x), "\\.")[[1]][1]}))
}

# extract all novel contigs
# as in the DE step, get all contigs that have
# an EC containing no reference transcripts
tx_ec <- data.table(distinct(ec_matrix[,c("ec_names", "transcript")]))
tx_ec$novel <- !tx_ec$transcript%in%txs
novel_contig_ecs <- tx_ec[, all(novel), by="ec_names"]
novel_contig_ecs <- unique(novel_contig_ecs$ec_names[novel_contig_ecs$V1])
novel_contigs <- unique(tx_ec$transcript[tx_ec$ec_names%in%novel_contig_ecs])

#############################################################
# Calculate VAFs
#############################################################

print("Estimating VAFs...")

# calculate the wildtype TPM by summing TPMs of all wildtype
# transcripts (anything that isn't a novel contig) per gene
qn <- data.frame(TPM=txi$abundance[,1], transcript=rownames(txi$abundance))
like_refseq <- qn$transcript%like%"hg38_ncbiRefSeq"
if(any(like_refseq)) {
    qn[like_refseq, 'transcript'] <- sapply(qn[like_refseq, 'transcript'],
                                            function(x){strsplit(gsub("hg38_ncbiRefSeq_", "", x), "\\.")[[1]][1]})
}
x <- inner_join(qn, tx2g, by="transcript")
wt_count <- data.table(x[!x$transcript%in%novel_contigs,])
wt_count <- distinct(wt_count)[, list(WT=sum(TPM, na.rm=TRUE)), by="gene"]
wt_count <- wt_count[wt_count$WT > 0,]

# now add the wildtype counts back to the quant table
# and extract only novel contigs
x <- left_join(x, wt_count, by="gene")
x <- x[x$transcript%in%cinfo$contig_id,]

# if contigs span multiple genes, we need to get the mean TPM
mean_tpm <- data.table(x)[, list(mean_WT_TPM=mean(WT, na.rm=TRUE)), by="transcript"]
x <- inner_join(x, mean_tpm, by=c("transcript"))
x$VAF <- x$TPM / (x$TPM + x$mean_WT_TPM)
x$VAF[is.nan(x$VAF)] <- 1 #assume no WT counts have a VAF of 1
colnames(x)[2] <- 'contig_id'
x <- x[,c('contig_id', 'gene', 'TPM', 'WT', 'mean_WT_TPM', 'VAF')]
if (nrow(x) > 0) {
    write.table(x, file=outfile, row.names=FALSE, quote=FALSE, sep="\t")
} else {
    stop("ERROR: no variants to output. Please double-check your reference files.")
}
