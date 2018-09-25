suppressWarnings(suppressMessages(require(edgeR)))
suppressMessages(require(reshape2))
suppressMessages(require(dplyr))
suppressMessages(require(data.table))
suppressMessages(require(IRanges))
suppressMessages(require(UpSetR))
suppressMessages(require(DEXSeq))
suppressMessages(require(stringr))
suppressMessages(require(seqinr))
options(stringsAsFactors = FALSE)

# load helper methods
args <- commandArgs(trailingOnly=FALSE)
file.arg <- grep("--file=",args,value=TRUE)
incl.path <- gsub("--file=(.*)compare_eq_classes.R","\\1helper_methods.R", file.arg)
source(incl.path, chdir=TRUE)

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 5) {
    args <- c("--help")
}

if("--help" %in% args) {
    cat("
        Load in Salmon equivalence classes from two samples,
        match them and perform differential splice analysis
        using equivalence classes as exons.

        Usage:
        Rscript compare_eq_classes.R <case_name> <ec_matrix> <annotation_file> <output> --iters=<n_iters> --sample=<N>\n\n

        Default iterations is 1 and default sample number is equivalent to number of controls
        Note: at least two control samples are required.\n\n")

    q(save="no")
}

FDR_cutoff <- 0.05
novel_contig_regex <- '^k[0-9]+_[0-9]+'

# parse arguments
args <- args[c(grep('--', args, invert=T), grep('--', args))]
case_name <- args[1]
ec_matrix_file <- args[2]
annotation_file <- args[3]
outfile <- args[4]

n_sample <- grep("--sample=",args, value=TRUE)
n_iters <- grep("--iters=",args, value=TRUE)
if (length(n_iters)==0) {
    n_iters <- 1
} else {
    n_iters <- strsplit(n_iters, "=")[[1]][2]
}
n_iters <- as.numeric(n_iters)

#############################################################
# load data
#############################################################
print('Reading input data...')
ec_matrix <- fread(ec_matrix_file, sep='\t')
annotation <- read.fasta(annotation_file, as.string=TRUE)

n_controls <- ncol(ec_matrix)-4 # all cols minus gene/tx/ec and case cols
if (length(n_sample)==0) {
    n_sample <- n_controls
} else {
    n_sample <- min(strsplit(n_sample, "=")[[1]][2], n_controls)
}
n_sample <- as.numeric(n_sample)

#############################################################
# annotation + prepare data for DE analysis
#############################################################
print('Generating transcript to gene lookup from fasta file...')
genes <- annotation %>% lapply(function(x){attributes(x)$Annot}) %>% str_match("gene_symbol:([a-zA-Z0-9.-]+)")
genes_tx <- distinct(data.frame(tx_id=names(annotation), gene=genes[,2]))
rm(genes, annotation); gc() #cleanup

print('Preparing expression table...')
# match ECs to genes
info <- match_tx_to_genes(ec_matrix, genes_tx)

#create reference lookup
tx_ec_gn <- info[,c('ec_names', 'genes', 'transcript')]

# get genes that are associated with all interesting novel contigs
tx_ec <- data.table(distinct(tx_ec_gn[,c('ec_names','transcript')]))
tx_to_ecs <- tx_ec[, paste(transcript, collapse=':'), by=list(ec_names)]
tx_to_ecs <- tx_to_ecs[grep('ENST', tx_to_ecs$V1, invert = T),] # ECs containing only novel contigs
colnames(tx_to_ecs)[2] <- 'contigs'
uniq_ecs <- tx_to_ecs$ec_names

## ambiguous mapping info, useful for final output
#ec_path <- paste(salmon_outdir, 'eq_classes.txt', sep='/')
#ambig_info_path <- paste(salmon_outdir, 'ambig_info.tsv', sep='/')
#uac <- get_ambig_info(ec_path, ambig_info_path, tx_ec_gn)

# cleanup
rm(ec_matrix); gc()

#############################################################
# perform DE
#############################################################

print('Performing differential expression analysis...')
bs_results <- run_edgeR(case_name, info, uniq_ecs, tx_to_ecs, dirname(outfile), cpm_cutoff=0)

#############################################################
# compile and write results
#############################################################

print('Compiling and writing results...')
concat_results <- bs_results
#concat_results <- left_join(concat_results, uac, by=c('ec_names','genes','transcript'))
#concat_results <- concat_results[order(concat_results$FDR),]

if (nrow(concat_results) > 0) {
    write.table(concat_results, outfile, row.names=F, quote=F, sep='\t')
} else {
    print('No significant ECs associated with novel contigs! Pipeline cannot continue.')
}
