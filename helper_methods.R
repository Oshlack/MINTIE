library(dplyr)

match_tx_to_genes <- function(ec_matrix, grp_novel, genes_tx) {
    genes_tx <- data.frame(genes_tx)
    genes_tx <- unique(genes_tx[,c('tx_id', 'symbol')])

    int_ec_matrix <- left_join(ec_matrix, grp_novel, by='transcript')
    int_ec_matrix <- within(int_ec_matrix, tx_id <- substr(transcript, 1, 15)) # remove dots in tx IDs

    int_ec_matrix <- left_join(int_ec_matrix, genes_tx, by='tx_id')
    int_ec_matrix[is.na(int_ec_matrix$gene), 'gene'] <- int_ec_matrix[is.na(int_ec_matrix$gene), 'symbol']

    ntx <- nrow(distinct(int_ec_matrix[,c('transcript', 'gene')]))
    ngs <- length(unique(int_ec_matrix$gene))
    print(paste(ntx, 'transcripts across', ngs, 'genes found'))

    # filter out any ECs with no associated gene
    int_ec_matrix <- int_ec_matrix[!is.na(int_ec_matrix$gene),]
    int_ec_matrix <- int_ec_matrix[int_ec_matrix$gene!='',]

    print('Grouping together all genes involved in fusions...')
    # collapse all genes involved in fusion into one 'super' gene
    fusion_contigs <- unique(int_ec_matrix[grep('\\|', int_ec_matrix$gene), 'gene'])
    fusion_contigs <- sapply(fusion_contigs, strsplit, split='\\|')
    rename_gns <- sapply(fusion_contigs, function(x){which(int_ec_matrix$gene%in%x)})
    for(i in 1:length(rename_gns)) {
        int_ec_matrix[rename_gns[[i]], 'gene'] <- names(rename_gns)[i]
    }

    ntx2 <- nrow(distinct(int_ec_matrix[,c('transcript', 'gene')]))
    info <- int_ec_matrix[, grep('cancer|control|ec|gene|transcript', colnames(int_ec_matrix))]

    print(paste('could not find genes for ', ntx-ntx2, ' transcripts (', round((ntx-ntx2)/ntx,4)*100 , '%)', sep=''))
    return(info)
}

get_de_novel_contigs <- function(top, tx_to_ecs, FDR_cutoff=0.05) {
    sig_ecs <- top[top$FDR<FDR_cutoff,]$ec_names
    x <- tx_to_ecs[tx_to_ecs$Group.1%in%sig_ecs,]
    x <- x[grep('ENST', x$transcript, invert = T), ]
    colnames(x)[1:2] <- c('ec_names', 'gene')
    x <- left_join(x, top, by=c('gene','ec_names'))
    x <- x[order(x$FDR, decreasing = F),]
    return(x)
}

is_interesting <- function(novel_contig, tx_blat) {
    # to be interesting, contig either contains multiple genes
    # or has gaps in all its reference transcriptome blast hits
    z <- data.table(tx_blat[tx_blat$V10==novel_contig,])
    if(length(unique(z$symbol))>1){return(TRUE)}
    z <- z[z$V8<=7,]
    if(nrow(z)==0){return(TRUE)}
    zr <- reduce(IRanges(z$V12, z$V13))
    cover_full_len <- z[(z$V13-z$V12)==z$V11]$V6
    gaps_in_reftx_whole_txs <- if(length(cover_full_len)>0){all(cover_full_len>7)}else{TRUE}
    gaps_in_reftx_concat_exons <- any(!zr@width-1 == unique(z$V11))
    return(gaps_in_reftx_concat_exons | gaps_in_reftx_whole_txs)
}

get_ambig_info <- function(ec_path, ambig_path, tx_ec_gn) {
    # get ambiguous mapping info from salmon,
    # match with transcripts and calculate
    # ambiguous mapping fraction
    x <- scan(ec_path, what = '', sep = '\n', quiet = T)
    x <- strsplit(x, '[ \t]+')
    max.col <- max(sapply(x, length))

    cn <- paste('X', 1:max.col, sep='')
    x <- read.delim(ec_path, header=F, sep='\t', col.names = cn)

    qs <- grep('^\\d+$',x$X1)[-c(1:2)] # get indexes of EC list
    txs <- x$X1[-c(1:2,qs)]

    y <- read.delim(ambig_path,sep='\t')
    uac <- cbind(transcript=txs, y)
    uac <- left_join(uac, tx_ec_gn, by='transcript')

    uac <- uac[(uac$AmbigCount+uac$UniqueCount)!=0,]
    uac$ambig_ratio <- uac$AmbigCount / (uac$AmbigCount +  uac$UniqueCount)

    return(uac)
}

bootstrap_diffsplice <- function(full_info, int_genes, n_controls, n_iters, select_ecs, tx_to_ecs) {
    # perform differential splicing selecting n_controls
    # bootstrap for n_iters iterations

    tx_ec_gn <- full_info[,c('ec_names', 'gene', 'transcript')] #create reference lookup
    info <- distinct(full_info[,grep('cancer|control|ec|gene', colnames(full_info))])
    info <- info[info$gene%in%int_genes,] #consider only ECs mapping to interesting genes

    print('Performing bootstrapped differential splicing analysis')
    results <- NULL
    for (i in 1:n_iters) {
        print(paste('Iteration', i, 'of', n_iters))

        # sample controls for comparison
        controls <- colnames(info)[grep('control', colnames(info))]
        if(n_controls<length(controls)){controls <- sample(controls, n_controls)}

        counts <- info[, c('cancer', controls)]
        counts <- as.matrix(apply(counts, 2, as.numeric))
        genes <- info[,c('gene', 'ec_names')]

        # CPM filtering
        keep <- rowSums(cpm(counts) > 1) >= 1
        counts <- counts[keep,]
        genes <- genes[keep,]

        # construct DGE object
        group <- factor(c('cancer', rep('control', ncol(counts)-1)))
        dge <- DGEList(counts=counts, group=group, genes=genes)
        dge <- calcNormFactors(dge)
        des <- model.matrix(~0+group)
        dge <- estimateDisp(dge, des)

        # fit and perform diff splice analysis
        fit <- glmFit(dge, des)
        sp <- diffSpliceDGE(fit, geneid='gene', exonid='ec_names', contrast=c(1,-1), verbose = FALSE)

        top_ds <- data.frame(topSpliceDGE(sp, test='Simes', n=Inf)) #gene-level diff splicing using Simes adjustment
        top_ex <- data.frame(topSpliceDGE(sp, test='exon', n=Inf)) #exon-level diff splicing

        sig_genes <- top_ds[top_ds$FDR<0.05, 'gene']
        sig_exons <- top_ex[top_ex$FDR<0.05, 'ec_names']

        # construct results dataframe
        txs_in_ec <- data.frame(table(tx_ec_gn$ec_names))
        txs_in_ec$Var1 <- as.character(txs_in_ec$Var1)
        colnames(txs_in_ec) <- c('ec_names', 'n_txs_in_ec')

        colnames(top_ds)[3:4] <- paste('Simes', colnames(top_ds)[3:4], sep='.')
        spg <- sp$genes[sp$genes$gene%in%sig_genes & sp$genes$ec_names%in%select_ecs,]
        spg <- left_join(spg, txs_in_ec, by='ec_names')
        spg <- left_join(spg, tx_to_ecs, by='ec_names')
        spg <- left_join(spg, top_ex, by=c('gene', 'ec_names'))
        spg <- left_join(spg, top_ds, by=c('gene'))

        spg$significant_ec <- spg$FDR<0.05
        nsig <- data.table(spg[,c('gene','significant_ec')])[,sum(significant_ec),by=gene]
        ntot <- data.table(spg[,c('gene','significant_ec')])[,length(significant_ec),by=gene]
        spg <- left_join(spg, nsig, by='gene')
        spg <- left_join(spg, ntot, by='gene')

        colnames(spg)[(ncol(spg)-1):ncol(spg)] <- c('sig_ecs_in_gene', 'total_ecs_in_gene')
        spg <- spg[,!colnames(spg)%in%'significant_ec']
        spg <- spg[spg$FDR<0.05 & spg$Simes.FDR<0.05,]

        results[[i]] <- list(controls=controls, spg=spg)
    }
    return(results)
}

concatenate_bs_results <- function(bs_results, n_iters) {
    # take results from bootstrap_diffsplice and return
    # gene/ECs present in all runs with mean FDRs and logFCs
    if(n_iters == 1) {
        concat_results <- bs_results[[1]]$spg
    } else {
        tmp <- NULL
        for(i in 1:n_iters) {
            tmp[[i]] <- bs_results[[i]]$spg
        }
        shared_cols <- c('gene','ec_names')
        ix_results <- Reduce(function(x, y) inner_join(x, y, by=shared_cols), tmp, accumulate = F)
        concat_results <- ix_results[,1:2]
        concat_results$n_txs_in_ec <- apply(ix_results[,grep('n_txs',colnames(ix_results))], 1, max)
        concat_results$NExons <- apply(ix_results[,grep('NExons',colnames(ix_results))], 1, max)
        concat_results$sig_ecs_in_gene <- apply(ix_results[,grep('sig_ecs_in_gene',colnames(ix_results))], 1, max)
        concat_results$total_ecs_in_gene <- apply(ix_results[,grep('total_ecs_in_gene',colnames(ix_results))], 1, max)
        concat_results$mean_FDR <- apply(ix_results[,grep('^FDR',colnames(ix_results))], 1, mean)
        concat_results$mean_Simes_FDR <- apply(ix_results[,grep('Simes.FDR',colnames(ix_results))], 1, mean)
        concat_results$mean_logFC <- sapply(apply(apply(ix_results[,grep('log',colnames(ix_results))],1,exp),2,mean),log)
        concat_results <- concat_results[order(concat_results$mean_FDR),]
    }
    return(concat_results)
}
