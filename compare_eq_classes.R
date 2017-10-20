suppressWarnings(suppressMessages(require(edgeR)))
suppressMessages(require(reshape2))

args <- commandArgs(TRUE)

if(length(args) < 3) {
    args <- c("--help")
}

if("--help" %in% args) {
    cat("
        Load in Salmon equivalence classes from two samples,
        match them and perform differential splice analysis
        using equivalence classes as exons.

        Usage:
        ./compare_eq_classes.R <ec_class1> <ec_class2> <output>\n\n")
    
    q(save="no")
}

# TODO: make these variables either auto-inferred or CMD-line arguments
bcv <- 0.1
FDR_cutoff <- 0.05

# load equivalence class file and return
# transcripts matched to equivalence classes
# with their respective counts
load_ecs <- function(fpath) {
    x <- scan(fpath, what = '', sep = '\n', quiet = T)
    x <- strsplit(x, '[ \t]+')
    max.col <- max(sapply(x, length))
    
    cn <- paste('X', 1:max.col, sep='')
    x <- read.delim(fpath, header=F, sep='\t', col.names = cn, stringsAsFactors = F)
    
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
    tx_ec <- data.frame(do.call(rbind, tx_ec), stringsAsFactors = F)
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
    xc <- data.frame(transcript=xc$transcript, xc_tmp, stringsAsFactors = F)
    rownames(xc) <- xc$transcript
    return(xc)
}

# create lookup table for matching equivalence classes
create_ec_lookup <- function(x, y) {
    xc <- create_ec_matrix(x)
    yc <- create_ec_matrix(y)
    common <- intersect(xc$transcript, yc$transcript)
    
    # create equivalence class transcript patterns
    xc_ecs <- as.character(apply(xc[common,2:ncol(xc)], 2, paste, collapse=''))
    yc_ecs <- as.character(apply(yc[common,2:ncol(yc)], 2, paste, collapse=''))
    
    ec_lookup <- cbind(ec_x=colnames(xc[common,2:ncol(xc)]), pattern=xc_ecs)
    ec_tmp <- cbind(ec_y=colnames(yc[common,2:ncol(yc)]), pattern=yc_ecs)
    ec_lookup <- merge(ec_lookup, ec_tmp, by='pattern')
    
    colnames(x)[2:3] <- c('ec_x', 'counts_x')
    colnames(y)[2:3] <- c('ec_y', 'counts_y')
    ec_lookup <- merge(x, ec_lookup, by='ec_x', all.x=T)
    ec_lookup <- merge(y, ec_lookup, by='ec_y', all.x=T)
    
    return(ec_lookup)
}
 
ec1_file <- args[1]
ec2_file <- args[2]

x <- load_ecs(ec1_file)
y <- load_ecs(ec2_file)
ec_lookup <- create_ec_lookup(x, y)

counts <- data.frame(x=ec_lookup$counts_x, y=ec_lookup$counts_y, stringsAsFactors = F)
counts[is.na(counts)] <- 0
counts <- as.matrix(apply(counts, 2, as.numeric))
genes <- data.frame(transcript=ec_lookup$transcript.x, ec=ec_lookup$ec_x, stringsAsFactors = F)
genes[is.na(genes$ec),'ec'] <- paste(ec_lookup[is.na(genes$ec),]$ec_y,'y',sep='_')

# CPM filtering
keep <- apply(cpm(counts), 1, sum) > 1
counts <- counts[keep,]
genes <- genes[keep,]

# make DGE object
group <- factor(c('x', 'y'))
dge <- DGEList(counts=counts, genes=genes, group=group)
dge <- calcNormFactors(dge)
des <- model.matrix(~ group)

# fit and perform differential splicing
fit <- glmFit(dge, des, robust=TRUE, dispersion = bcv^2)
sp <- diffSpliceDGE(fit, geneid='transcript', exonid='ec', verbose = F)
top <- topSpliceDGE(sp, test="Simes", n=Inf)

outfile <- args[3]
output <- top[top$FDR<FDR_cutoff,]
write.table(output, outfile, row.names=F, quote=F, sep='\t')
