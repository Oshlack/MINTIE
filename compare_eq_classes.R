suppressMessages(require(edgeR))
suppressMessages(require(reshape2))
suppressMessages(require(dplyr))
options(stringsAsFactors = FALSE)

args <- commandArgs(TRUE)

if(length(args) < 4) {
    args <- c("--help")
}

if("--help" %in% args) {
    cat("
        Load in Salmon equivalence classes from two samples,
        match them and perform differential splice analysis
        using equivalence classes as exons.

        Usage:
        ./compare_eq_classes.R <cancer_ec_class> <control1_ec> <control2_ec> ... <output>\n\n

        Note: at least two control samples are required.\n\n")

    q(save="no")
}

# TODO: make these variables either auto-inferred or CMD-line arguments
FDR_cutoff <- 0.05

# load equivalence class file and return
# transcripts matched to equivalence classes
# with their respective counts
load_ecs <- function(fpath) {
    x <- scan(fpath, what = '', sep = '\n', quiet = T)
    x <- strsplit(x, '[ \t]+')
    max.col <- max(sapply(x, length))

    cn <- paste('X', 1:max.col, sep='')
    x <- read.delim(fpath, header=F, sep='\t', col.names = cn)

    qs <- grep('^\\d+$',x$X1)[-c(1:2)] # get indexes of EC list
    txs <- x$X1[-c(1:2,qs)]

    y <- x[qs,]
    tx_ec <- NULL
    for(i in 1:nrow(y)) { # add trancript/ec associations
        ntxs <- as.numeric(y[i,1]) # number of TXs in this EC
        ec <- paste('ec', i, sep='')
        for(j in 1:ntxs) {
            tx_id <- y[i, j+1] + 1
            tx <- txs[tx_id]
            count <- y[i, ntxs+2]
            tx_ec[[paste(tx,ec)]] <- c(tx, ec, count)
        }
    }
    tx_ec <- data.frame(do.call(rbind, tx_ec))
    colnames(tx_ec) <- c('transcript', 'ec', 'count')
    rownames(tx_ec) <- 1:nrow(tx_ec)

    return(tx_ec)
}

# transform equivalence class list into
# transcript presence/absence matrix
create_ec_matrix <- function(x) {
    xc <- dcast(x[,1:2], transcript ~ ec, value.var = 'ec')
    xc_tmp <- xc[,2:ncol(xc)]
    xc_tmp[!is.na(xc_tmp)] <- 1
    xc_tmp[is.na(xc_tmp)] <- 0

    xc_tmp <- apply(xc_tmp, 2, as.numeric)
    xc <- data.frame(transcript=xc$transcript, xc_tmp)
    rownames(xc) <- xc$transcript
    return(xc)
}

# create lookup table for matching equivalence classes
create_ec_lookup <- function(x, all_tx) {
    xc <- create_ec_matrix(x)
    xc <- xc[all_tx,]; xc[is.na(xc)] <- 0
    xc_ecs <- as.character(apply(xc[,2:ncol(xc)], 2, paste, collapse=''))
    ec_lookup <- data.frame(pattern=xc_ecs, ec=colnames(xc[,2:ncol(xc)]))
    ec_lookup <- inner_join(ec_lookup, x, by='ec')
    return(ec_lookup)
}

ec_cancer_file <- args[1]
control_files <- args[2:(length(args)-1)]

x <- load_ecs(ec_cancer_file)

all_tx <- list()
all_tx[[1]] <- x$transcript
controls <- list()
for(i in 1:length(control_files)) {
    control_file <- control_files[i]
    controls[[i]] <- load_ecs(control_file)
    all_tx[[i+1]] <- unique(controls[[i]]$transcript)
}
all_tx <- Reduce(union, all_tx)

lookup <- list()
lookup[[1]] <- create_ec_lookup(x, all_tx)
for(i in 1:length(controls)) {
    con <- create_ec_lookup(controls[[i]], all_tx)
    con <- con[,c('pattern', 'transcript', 'count')]
    colnames(con)[3] <- paste('control',i,sep='')
    lookup[[i+1]] <- con
}

ec_union <- Reduce(function(x, y) full_join(x, y, by=c('pattern','transcript')), lookup, accumulate = F)
info <- unique(ec_union)
genes <- ec_union[,c('transcript', 'ec')]

counts <- info[,-c(1:3)]
counts[is.na(counts)] <- 0
counts <- apply(counts, 2, as.numeric)

#make new equivalence class labels
max_ec <- max(as.numeric(sapply(unique(genes$ec), function(x){strsplit(x,"ec")[[1]][2]})), na.rm=T)
new_ec <- paste('ec', (max_ec+1):(sum(is.na(genes$ec))+max_ec), sep='')
genes[is.na(genes$ec),'ec'] <- new_ec

# CPM filtering
keep <- rowSums(cpm(counts) > 1) >= ncol(counts)/3
keep <- apply(cbind(!is.na(genes$transcript), keep), 1, all)
counts <- counts[keep,]
genes <- genes[keep,]

# make DGE object, estimate dispersion
group <- factor(c('cancer', rep('control', ncol(counts)-1)))
dge <- DGEList(counts=counts, genes=genes, group=group)
dge <- calcNormFactors(dge)
des <- model.matrix(~ group)
dge <- estimateDisp(dge, des, robust=TRUE)

# fit and perform differential splicing
fit <- glmQLFit(dge, des, robust=TRUE)
qlf <- glmQLFTest(fit)
sp <- diffSpliceDGE(fit, geneid='transcript', exonid='ec', verbose = F, contrast=c(-1,1))
top <- topSpliceDGE(sp, test='Simes', n=Inf)

# filter and write output
outfile <- args[length(args)]
output <- top[top$FDR<FDR_cutoff,]
write.table(output, outfile, row.names=F, quote=F, sep='\t')
