library(dplyr)
library(data.table)
library(edgeR)

run_edgeR <- function(case_name, ec_matrix, tx_ec, outdir, cpm_cutoff=0.1, qval=0.05, min_logfc=0) {
    data <- distinct(ec_matrix[, !colnames(ec_matrix)%in%c("transcript", "tx_ids"), with=F])
    data <- data[data$ec_names%in%tx_ec$ec_names,]
    if(nrow(data) == 0) {
        stop("No novel ECs exist in EC matrix. Please check EC matrix input.")
    }

    # prepare counts matrix
    counts <- data[,!colnames(data)%in%"ec_names", with=F]
    counts <- as.matrix(apply(counts, 2, as.numeric))
    rownames(counts) <- data$ec_names

    # define groups
    group <- rep("control", ncol(counts))
    group[colnames(counts)==case_name] <- "case"
    if(!"case"%in%group){
        group[colnames(counts)==gsub("-", ".", case_name)] <- "case"
    }
    if(!"case"%in%group) {
        group[colnames(counts)==paste0("X", case_name)] <- "case"
    }
    if(!"case"%in%group) {
        group[colnames(counts)==gsub("-", ".", paste0("X", case_name))] <- "case"
    }
    if(!"case"%in%group) {
        print(paste("Case name:", case_name))
        print(paste("All sample names:", paste(colnames(counts), collapse=", ")))
        stop("No case sample found. Please ensure case name matches sample columns.")
    }
    group <- factor(group, levels=c("control", "case"))
    colnames(counts) <- as.character(sapply(colnames(counts), function(x){gsub("-", "_", x)}))
    des <- model.matrix(~group)
    colnames(des)[2] <- "case"

    # filter
    keep <- as.numeric(cpm(counts)[, group=="case"]) > cpm_cutoff
    counts <- counts[keep,]
    if(nrow(counts) == 0) {
        stop(paste("All ECCs were below the CPM cutoff of", cpm_cutoff, "."))
    }

    # make counts summary
    case_counts <- counts[, group=="case"]
    con_counts <- data.frame(counts[, group=="control"])
    counts_summary <- data.frame(ec_names=names(case_counts),
                                 case_reads=case_counts,
                                 controls_total_reads=apply(con_counts, 1, sum))

    # write counts table to file for reference
    counts_out <- data.frame(ec_names=rownames(counts), counts)
    write.table(counts_out, file=paste(outdir, "counts_table.txt", sep="/"), row.names=F, quote=F, sep="\t")

    # perform DE
    dge <- DGEList(counts = counts, group = group)
    dge <- calcNormFactors(dge)
    dge <- estimateGLMCommonDisp(dge, design=des, verbose=T)
    dge <- estimateGLMTrendedDisp(dge, design=des)
    dge <- estimateGLMTagwiseDisp(dge, design=des)

    fit <- glmQLFit(dge, design=des)
    qlf <- glmQLFTest(fit, coef=2)
    dx_df <- data.frame(topTags(qlf, n=Inf))
    significant <- dx_df$FDR<qval & dx_df$logFC>min_logfc
    if(nrow(dx_df[significant,]) == 0) {
        stop("No significantly differentially expressed ECCs were found. Pipeline cannot continue.")
    }
    print("Successfully ran differential expression!")
    print(paste("There were", sum(significant), "significantly differentially expressed novel-contig ECs."))

    # get number of contigs in each EC
    contigs_in_ec <- data.frame(table(tx_ec$ec_names))
    contigs_in_ec$Var1 <- as.character(contigs_in_ec$Var1)
    colnames(contigs_in_ec) <- c("ec_names", "n_contigs_in_ec")

    dx_df$ec_names <- rownames(dx_df)
    dx_df <- left_join(dx_df, contigs_in_ec, by="ec_names")
    dx_df <- left_join(dx_df, tx_ec, by="ec_names")
    dx_df <- left_join(dx_df, counts_summary, by="ec_names")
    dx_df <- dx_df[order(dx_df$PValue, decreasing=F),]

    # write full results
    write.table(dx_df, file=paste(outdir, "full_edgeR_results.txt", sep="/"), row.names=F, quote=F, sep="\t")
    dx_df <- dx_df[dx_df$FDR<qval & dx_df$logFC>min_logfc,]

    return(dx_df)
}
