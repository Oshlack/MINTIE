suppressMessages(require(edgeR))
suppressMessages(require(reshape2))
suppressMessages(require(dplyr))
options(stringsAsFactors = FALSE)

args <- commandArgs(TRUE)

if(length(args) < 3) {
    args <- c("--help")
}

if("--help" %in% args) {
    cat("
        Perform differential expression and
        differential splice analysis using
        equivalence classes.

        Usage:
        Rscript compare_eq_classes.R <ec_matrix> <groupings> <output>\n\n

        Note: at least two control samples are required.\n\n")

    q(save="no")
}

# TODO: make these variables either auto-inferred or CMD-line arguments
FDR_cutoff <- 0.05
novel_contig_regex <- '^k[0-9]+_[0-9]+'

ec_matrix_file <- args[1]
groupings_file <- args[2]
output_prefix <- args[3]

ec_matrix <- read.delim(ec_matrix_file, sep='\t')

################## diffsplice testing using ECs ##################

# parse contig to gene groupings
grp <- read.delim(groupings_file, sep='\t', header=F)
colnames(grp) <- c('transcript','contig')
grp_novel <- grp[grep(novel_contig_regex, grp$transcript),]
grp_gene <- grp[grep(novel_contig_regex, grp$transcript, invert = T),]

# create gene-contig lookup
contig_gene <- merge(grp_novel, grp_gene, by='contig')[,-1]
colnames(contig_gene) <- c('transcript', 'gene')

# remove denovo assembled contigs that don't group with any genes
tmp <- merge(ec_matrix, contig_gene, all.x=T, by='transcript')
tmp[is.na(tmp$gene),'gene'] <- tmp[is.na(tmp$gene),'transcript']
tmp <- tmp[grep(novel_contig_regex, tmp$gene, invert = T),]
info <- unique(tmp[, grep('cancer|control|ec|gene', colnames(tmp))])

# make counts and genes dataframes
counts <- info[, grep('cancer|control', colnames(info))]
counts <- as.matrix(apply(counts, 2, as.numeric))
genes <- info[,c('gene', 'ec_names')]

# filter low count ECs
keep <- rowSums(cpm(counts) > 1) >= 1
counts <- counts[keep,]
genes <- genes[keep,]

# make DGE object, estimate dispersion
group <- factor(c('cancer', rep('control', ncol(counts)-1)))
dge <- DGEList(counts=counts, group=group, genes=genes)
dge <- calcNormFactors(dge)
des <- model.matrix(~ group)
dge <- estimateDisp(dge, des, robust=TRUE)

# perform diffsplice on equivalance classes of genes
fit <- glmFit(dge, des)
sp <- diffSpliceDGE(fit, contrast = c(1,-1), geneid='gene', exonid='ec_names')
top <- topSpliceDGE(sp, n=Inf)
top <- top[top$FDR<FDR_cutoff,]

# write output
outfile <- paste(output_prefix, 'diffsplice.txt', sep='_')
write.table(top, outfile, row.names=F, quote=F, sep='\t')

################## DE using equivalence classes DE ##################

# keep interesting ECs containing de novo assemblies
interesting_ecs <- unique(ec_matrix[grep(novel_contig_regex, ec_matrix$transcript), 'ec_names'])
ec_matrix <- ec_matrix[ec_matrix$ec_names%in%interesting_ecs,]

# make counts
counts <- unique(ec_matrix[, grep('ec|cancer|control', colnames(ec_matrix))])
ecs <- counts$ec_names
counts <- apply(counts[,-grep('ec', colnames(ec_matrix))], 2, as.numeric)
rownames(counts) <- ecs

# filter low count ECs
keep <- rowSums(cpm(counts) > 1) >= 1
counts <- counts[keep,]

# make DGE object, estimate dispersion
group <- factor(c('cancer', rep('control', ncol(counts)-1)))
dge <- DGEList(counts=counts, group=group)
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge, des, robust=TRUE)

# fit and perform differential splicing
fit <- glmQLFit(dge, des, robust=TRUE)
qlf <- glmQLFTest(fit, contrast=c(1,-1))
top <- data.frame(topTags(qlf, n=Inf))
top <- data.frame(ec_names=rownames(top), top)
top <- top[top$FDR<FDR_cutoff,]

# filter and write output
ec_matrix_out <- ec_matrix[ec_matrix$ec_names%in%top$ec_names,grep('cancer|control|ec', colnames(ec_matrix))]
output <- left_join(ec_matrix_out, top, by='ec_names')
outfile <- paste(output_prefix, 'de.txt', sep='_')
write.table(output, outfile, row.names=F, quote=F, sep='\t')
