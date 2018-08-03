library(dplyr)
library(data.table)
library(edgeR)

match_tx_to_genes <- function(ec_matrix, grp_novel, genes_tx) {
    ec_matrix$tx_id <- ec_matrix$transcript
    int_ec_matrix <- left_join(ec_matrix, grp_novel, by='transcript')
    int_ec_matrix <- left_join(int_ec_matrix, genes_tx, by='tx_id')
    int_ec_matrix[is.na(int_ec_matrix$gene), 'gene'] <- int_ec_matrix[is.na(int_ec_matrix$gene), 'symbol']

    tmp <- distinct(int_ec_matrix[,c('transcript', 'gene')])
    ngs <- length(unique(int_ec_matrix$gene))
    print(paste(nrow(tmp), 'transcripts across', ngs, 'genes found'))

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

    ntx <- nrow(distinct(int_ec_matrix[,c('transcript', 'gene')]))
    info <- int_ec_matrix[, !colnames(int_ec_matrix)%in%c('tx_id', 'symbol')]

    if(nrow(tmp)-ntx!=0) {
        print(paste('could not find genes for ', nrow(tmp)-ntx, ' transcripts (',
                    round((nrow(tmp)-ntx)/nrow(tmp),4)*100 , '%):', sep=''))
        print(tmp$transcript[!tmp$transcript%in%info$transcript])
    }
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

#TODO: write method for extracting counts matrix
run_edgeR <- function(case_name, full_info, int_genes, select_ecs, tx_to_ecs, outdir, threads=8, cpm_cutoff=2, qval=0.05) {
    tx_ec_gn <- full_info[,c('ec_names', 'gene', 'transcript')] #create reference lookup
    info <- distinct(full_info[full_info$gene%in%int_genes, !colnames(full_info)%in%'transcript'])

    # collapse reference equivalence classes
    print('Collapsing reference equivalence classes')
    info$ec <- 'reference'
    info$ec[info$ec_names%in%select_ecs] <- info$ec_names[info$ec_names%in%select_ecs]
    info <- data.table(info[,!colnames(info)%in%'ec_names'])[, lapply(.SD, sum), by=c('gene','ec')]
    colnames(info)[colnames(info)=='ec'] <- 'ec_names'

    results <- NULL
    print('Performing differential expression analysis with edgeR...')
    start.time <- Sys.time(); print(start.time)

    counts <- data.frame(info)[,!colnames(info)%in%c('ec_names','gene')]
    counts <- as.matrix(apply(counts, 2, as.numeric))
    genes <- data.frame(info)[,c('gene', 'ec_names')]

    keep <- rowSums(cpm(counts) > cpm_cutoff) >= 1
    counts <- counts[keep,]
    genes <- genes[keep,]

    group <- rep('control', ncol(counts))
    group[colnames(counts)==case_name] <- 'cancer'
    if(!'cancer'%in%group){
        group[colnames(counts)==gsub('-','.',case_name)] <- 'cancer'
    }
    group <- factor(group, levels=c('control', 'cancer'))
    colnames(counts) <- as.character(sapply(colnames(counts), function(x){gsub('-', '_',x)}))

    dge <- DGEList(counts = counts, genes = paste(genes$ec_names, genes$gene), group = group)
    dge <- calcNormFactors(dge)

    des <- model.matrix(~group)
    colnames(des)[2] <- 'cancer'

    dge <- estimateGLMCommonDisp(dge, design=des, verbose=T)
    dge <- estimateGLMTrendedDisp(dge, design=des)
    dge <- estimateGLMTagwiseDisp(dge, design=des)
    disp <- dge$common.dispersion

    fit <- glmQLFit(dge, design=des)
    qlf <- glmQLFTest(fit, coef=2)

    top <- data.frame(topTags(qlf, n=Inf))
    tmp <- data.frame(t(data.frame(strsplit(top$genes, ' '))))
    colnames(tmp) <- c('ec_names', 'gene')
    dx_df <- cbind(top, tmp)
    dx_df <- merge(dx_df, tx_ec_gn, by=c('ec_names', 'gene'), all.x=T)

    txs_in_ec <- data.frame(table(tx_ec_gn$ec_names))
    txs_in_ec$Var1 <- as.character(txs_in_ec$Var1)
    colnames(txs_in_ec) <- c('ec_names', 'n_txs_in_ec')

    dx_df <- left_join(dx_df, txs_in_ec, by='ec_names')
    dx_df <- left_join(dx_df, tx_to_ecs, by='ec_names')

    # write full results
    write.table(dx_df, file=paste(outdir, 'full_edgeR_results.txt', sep='/'), row.names=F, quote=F, sep='\t')
    dx_df <- dx_df[dx_df$FDR<qval,]

    return(dx_df)
}

run_dexseq <- function(case_name, full_info, int_genes, select_ecs, tx_to_ecs, outdir, threads=8, cpm_cutoff=2) {
    tx_ec_gn <- full_info[,c('ec_names', 'gene', 'transcript')] #create reference lookup
    info <- distinct(full_info[full_info$gene%in%int_genes, !colnames(full_info)%in%'transcript'])
    print(head(info))

    # collapse reference equivalence classes
    print('Collapsing reference equivalence classes')
    info$ec <- 'reference'
    info$ec[info$ec_names%in%select_ecs] <- info$ec_names[info$ec_names%in%select_ecs]
    info <- data.table(info[,!colnames(info)%in%'ec_names'])[, lapply(.SD, sum), by=c('gene','ec')]
    colnames(info)[colnames(info)=='ec'] <- 'ec_names'

    results <- NULL
    print('Performing differential splicing analysis with DEXSeq...')
    start.time <- Sys.time(); print(start.time)

    counts <- data.frame(info)[,!colnames(info)%in%c('ec_names','gene')]
    counts <- as.matrix(apply(counts, 2, as.numeric))
    genes <- data.frame(info)[,c('gene', 'ec_names')]

    keep <- rowSums(cpm(counts) > cpm_cutoff) >= 1
    counts <- counts[keep,]
    genes <- genes[keep,]

    group <- rep('control', ncol(counts))
    #case_name <- gsub('-', '.', case_name)
    group[colnames(counts)==case_name] <- 'cancer'
    colnames(counts) <- as.character(sapply(colnames(counts), function(x){gsub('-', '_',x)}))

    sampleTable <- data.frame(condition=factor(group, levels=c('control', 'cancer')))
    rownames(sampleTable) <- colnames(counts)

    BPPARAM=MulticoreParam(workers=threads)
    dx <- DEXSeqDataSet(counts, sampleData = sampleTable,
                        groupID = genes$gene, featureID = genes$ec_names)
    dx <- estimateSizeFactors(dx)

    #print('Estimating dispersions...')
    #dx <- estimateDispersions(dx, BPPARAM=BPPARAM)
    #dx <- testForDEU(dx)
    #dx <- estimateExonFoldChanges(dx)

    print('Performing DEU test...')
    dxr <- DEXSeq(dx, BPPARAM=BPPARAM)
    pgq <- perGeneQValue(dxr)
    pgq <- data.frame(gene=names(pgq), gene.FDR=pgq)

    dx_df <- data.frame(dxr)
    print(head(dx_df))
    colnames(dx_df)[1:2] <- c('gene', 'ec_names')

    txs_in_ec <- data.frame(table(tx_ec_gn$ec_names))
    txs_in_ec$Var1 <- as.character(txs_in_ec$Var1)
    colnames(txs_in_ec) <- c('ec_names', 'n_txs_in_ec')

    dx_df <- left_join(dx_df, txs_in_ec, by='ec_names')
    dx_df <- left_join(dx_df, tx_to_ecs, by='ec_names')
    dx_df <- left_join(dx_df, pgq, by='gene')
    dx_df <- dx_df[, !colnames(dx_df)%in%'genomicData']

    dx_df$significant_ec <- dx_df$padj<0.05
    dx_df$significant_gene <- dx_df$gene.FDR<0.05
    nsig <- data.table(dx_df[,c('gene','significant_ec')])[,sum(significant_ec),by=gene]
    ntot <- data.table(dx_df[,c('gene','significant_ec')])[,length(significant_ec),by=gene]
    dx_df <- left_join(dx_df, nsig, by='gene')
    dx_df <- left_join(dx_df, ntot, by='gene')
    colnames(dx_df)[(ncol(dx_df)-1):ncol(dx_df)] <- c('sig_ecs_in_gene', 'total_ecs_in_gene')
    dx_df$novel_ec <- dx_df$ec_names%in%select_ecs

    end.time <- Sys.time(); print(end.time)
    time.taken <- end.time - start.time

    print(paste('Testing with DEXSeq took', time.taken, 'minutes!'))

    # write full results
    write.table(dx_df, file=paste(outdir, 'full_dexseq_results.txt', sep='/'), row.names=F, quote=F, sep='\t')

    #dx_df <- dx_df[!is.na(dx_df$significant_gene) & !is.na(dx_df$contigs) & !is.na(dx_df$significant_ec),]
    #dx_df <- dx_df[dx_df$significant_gene & dx_df$significant_ec & dx_df$novel_ec,]
    dx_df <- dx_df[!is.na(dx_df$contigs) & !is.na(dx_df$significant_ec),]
    dx_df <- dx_df[dx_df$significant_ec & dx_df$novel_ec,]

    return(dx_df)
}

bootstrap_diffsplice <- function(case_name, full_info, int_genes, n_controls, n_iters, select_ecs, tx_to_ecs, outdir, cpm_cutoff=2) {
    # perform differential splicing selecting n_controls
    # bootstrap for n_iters iterations
    tx_ec_gn <- full_info[,c('ec_names', 'gene', 'transcript')] #create reference lookup
    info <- distinct(full_info[,!colnames(full_info)%in%'transcript'])
    info <- info[info$gene%in%int_genes,] #consider only ECs mapping to interesting genes

    results <- NULL
    print('Performing differential splicing analysis...')
    for (i in 1:n_iters) {
        if(n_iters>1) {
            print(paste('Iteration', i, 'of', n_iters))
        }
        counts <- info[,!colnames(info)%in%c('ec_names','gene')]
        counts <- as.matrix(apply(counts, 2, as.numeric))
        genes <- info[,c('gene', 'ec_names')]

        # sample controls for comparison
        case_name <- gsub('-', '.', case_name)
        controls <- colnames(counts)[grep(case_name, colnames(counts), invert=T)]
        if(n_controls < length(controls)) {
            controls <- sample(controls, n_controls)
            counts <- counts[,c(case_name, controls)]
        }

        # CPM filtering
        keep <- rowSums(cpm(counts) > cpm_cutoff) >= 1
        counts <- counts[keep,]
        genes <- genes[keep,]

        # construct DGE object
        group <- rep('control', ncol(counts))
        group[colnames(counts)==case_name] <- 'cancer'
        colnames(counts) <- as.character(sapply(colnames(counts), function(x){gsub('-', '_',x)}))

        dge <- DGEList(counts=counts, group=as.factor(group), genes=genes)
        dge <- calcNormFactors(dge)
        des <- model.matrix(~0+as.factor(group))
        dge <- estimateDisp(dge, des)

        # fit and perform diff splice analysis
        fit <- glmFit(dge, des)
        sp <- diffSpliceDGE(fit, geneid='gene', exonid='ec_names', contrast=c(1,-1), verbose = FALSE)

        top_ds <- data.frame(topSpliceDGE(sp, test='gene', n=Inf)) #gene-level diff splicing
        top_ex <- data.frame(topSpliceDGE(sp, test='exon', n=Inf)) #exon-level diff splicing

        #sig_genes <- top_ds[top_ds$FDR<0.05, 'gene']
        #sig_exons <- top_ex[top_ex$FDR<0.05, 'ec_names']

        # construct results dataframe
        txs_in_ec <- data.frame(table(tx_ec_gn$ec_names))
        txs_in_ec$Var1 <- as.character(txs_in_ec$Var1)
        colnames(txs_in_ec) <- c('ec_names', 'n_txs_in_ec')

        colnames(top_ds)[4:5] <- paste('gene', colnames(top_ds)[4:5], sep='.')
        spg <- sp$genes
        spg <- left_join(spg, txs_in_ec, by='ec_names')
        spg <- left_join(spg, tx_to_ecs, by='ec_names')
        spg <- left_join(spg, top_ex, by=c('gene', 'ec_names'))
        spg <- left_join(spg, top_ds, by=c('gene'))

        spg$significant_ec <- spg$FDR<0.05
        spg$significant_gene <- spg$gene.FDR<0.05
        spg$novel_ec <- spg$ec_names%in%select_ecs

        nsig <- data.table(spg[,c('gene','significant_ec')])[,sum(significant_ec),by=gene]
        ntot <- data.table(spg[,c('gene','significant_ec')])[,length(significant_ec),by=gene]
        spg <- left_join(spg, nsig, by='gene')
        spg <- left_join(spg, ntot, by='gene')
        colnames(spg)[(ncol(spg)-1):ncol(spg)] <- c('sig_ecs_in_gene', 'total_ecs_in_gene')

        # write full results
        write.table(spg, file=paste(outdir, '/full_diffsplice_results_iter_', i, '.txt', sep=''), row.names=F, quote=F, sep='\t')

        spg <- spg[spg$significant_ec & spg$significant_gene & spg$novel_ec,]
        spg <- spg[,!colnames(spg)%in%c('significant_ec','significant_gene')]

        results[[i]] <- list(controls=controls, spg=spg)
    }
    return(results)
}

concatenate_bs_results <- function(bs_results, n_iters) {
    # take results from bootstrap_diffsplice and return
    # gene/ECs present in all runs with mean FDRs and logFCs
    if(n_iters == 1) {
        concat_results <- bs_results[[1]]$spg
        concat_results <- concat_results[order(concat_results$gene.FDR),]
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
        concat_results$mean_gene_FDR <- apply(ix_results[,grep('gene.FDR',colnames(ix_results))], 1, mean)
        concat_results$mean_logFC <- sapply(apply(apply(ix_results[,grep('log',colnames(ix_results))],1,exp),2,mean),log)
        concat_results <- concat_results[order(concat_results$mean_FDR),]
    }
    return(concat_results)
}
