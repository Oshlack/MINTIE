suppressWarnings(suppressMessages(require(edgeR)))
suppressMessages(require(reshape2))
suppressMessages(require(dplyr))
suppressMessages(require(EnsDb.Hsapiens.v86))
suppressMessages(require(data.table))
suppressMessages(require(IRanges))
suppressMessages(require(UpSetR))
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
        Rscript compare_eq_classes.R <ec_matrix> <groupings> <salmon_outdir> <tx_blat> <output> --iters=<n_iters> --sample=<N>\n\n

        Default iterations is 1 and default sample number is equivalent to number of controls
        Note: at least two control samples are required.\n\n")

    q(save="no")
}

FDR_cutoff <- 0.05
novel_contig_regex <- '^k[0-9]+_[0-9]+'

# parse arguments
args <- args[c(grep('--', args, invert=T), grep('--', args))]
ec_matrix_file <- args[1]
groupings_file <- args[2]
salmon_outdir <- args[3]
against_txome_blat <- args[4]
outfile <- args[5]

n_sample <- grep("--sample=",args,value=TRUE)
n_iters <- grep("--iters=",args,value=TRUE)
if (length(n_iters)==0) {
    n_iters <- 1
} else {
    n_iters <- strsplit(n_iters, "=")[[1]][2]
}
n_iters <- as.numeric(n_iters)

################## load data ##################
ec_matrix <- read.delim(gzfile(ec_matrix_file), sep='\t')
all.groupings <- read.delim(groupings_file, sep='\t', header=F)
colnames(all.groupings) <- c('transcript', 'gene')

n_controls <- length(grep('control', colnames(ec_matrix)))
if (length(n_sample)==0) {
    n_sample <- n_controls
} else {
    n_sample <- min(strsplit(n_sample, "=")[[1]][2], n_controls)
}
n_sample <- as.numeric(n_sample)

# match ECs to genes
grp_novel <- all.groupings[grep(novel_contig_regex, all.groupings$transcript),]
genes_tx <- transcripts(EnsDb.Hsapiens.v86, columns=listColumns(EnsDb.Hsapiens.v86, c('tx', 'gene')))
info <- match_tx_to_genes(ec_matrix, grp_novel, genes_tx)

print('Identifying novel contigs...')
tx_blat <- read.delim(against_txome_blat, sep='\t', skip = 5, header=F)
tx_blat$tx_id <- as.character(sapply(tx_blat$V14, function(x){strsplit(x, '\\.')[[1]][1]}))
tx_blat <- left_join(tx_blat, data.frame(genes_tx)[,c('tx_id', 'symbol')], by='tx_id')

tx_ec_gn <- info[,c('ec_names', 'gene', 'transcript')] #create reference lookup
nc <- unique(tx_ec_gn[grep(novel_contig_regex, tx_ec_gn$transcript),]$transcript) # all novel contigs
int_contigs <- sapply(nc, is_interesting, tx_blat) # novel contigs with consistent reference gaps
int_contigs <- names(int_contigs[as.logical(int_contigs)])

# get genes that are associated with all interesting novel contigs
int_ecs <- unique(ec_matrix[ec_matrix$transcript%in%int_contigs,]$ec_names)
tx_ec <- data.table(distinct(tx_ec_gn[,c('ec_names','transcript')]))
tx_to_ecs <- tx_ec[, paste(transcript, collapse=':'), by=list(ec_names)]
tx_to_ecs <- tx_to_ecs[grep('ENST',tx_to_ecs$V1, invert = T),] # ECs containing only novel contigs
colnames(tx_to_ecs)[2] <- 'contigs'
uniq_ecs <- tx_to_ecs$ec_names

# get interesting genes (any genes associated with non-reference ECs that only map to denovo contigs)
int_ecs <- intersect(uniq_ecs, int_ecs)
int_genes <- unique(info[info$ec_names%in%int_ecs,]$gene)

# ambiguous mapping info, useful for final output
ec_path <- paste(salmon_outdir, 'eq_classes.txt', sep='/')
ambig_info_path <- paste(salmon_outdir, 'ambig_info.tsv', sep='/')
uac <- get_ambig_info(ec_path, ambig_info_path, tx_ec_gn)

################## diffsplice testing using ECs ##################

bs_results <- bootstrap_diffsplice(info, int_genes, n_sample, n_iters, uniq_ecs, tx_to_ecs)

################## compile and write results ##################

print('Compiling and writing results...')
bs_genes <- NULL
for(i in 1:n_iters) {
    bs_genes[[i]] <- unique(bs_results[[i]]$spg$gene)
}

names(bs_genes) <- 1:n_iters
bs_up <- fromList(bs_genes)
if (n_iters > 1) {
    plotfile <- paste(dirname(outfile), '/upset_plot_', n_iters, '_iters.pdf', sep='')
    pdf(plotfile, width=8, height=5)
    upset(bs_up, order.by='freq', nsets=n_iters)
    dev.off()
}

concat_results <- concatenate_bs_results(bs_results, n_iters)
concat_results <- left_join(concat_results, uac, by=c('ec_names','gene'))
write.table(concat_results, outfile, row.names=F, quote=F, sep='\t')
