library(dplyr)
library(data.table)
library(edgeR)

run_edgeR <- function(case_name, ec_matrix, tx_ec, outdir, cpm_cutoff=0.1, qval=0.05, min_logfc=5) {
    data <- distinct(ec_matrix[, !colnames(ec_matrix)%in%c("transcript", "tx_ids", "tx_id"), with=FALSE])

    # prepare counts matrix
    counts <- data[,!colnames(data)%in%"ec_names", with=FALSE]
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

    # calculate "true" library sizes
    dge <- DGEList(counts = counts, group = group)
    dge <- calcNormFactors(dge)
    norm_factors <- dge$samples

    # filter
    case_cpms <- cpm(counts)[, group=="case"] %>% as.numeric()
    keep <- case_cpms > cpm_cutoff
    unfilt_counts <- counts
    counts <- counts[keep,]
    if(nrow(counts) == 0) {
        stop(paste("All ECCs were below the CPM cutoff of", paste0(cpm_cutoff, ".")))
    }
    counts <- counts[rownames(counts)%in%tx_ec$ec_names,]
    if(nrow(data) == 0) {
        stop("No novel ECs exist in EC matrix. Please check EC matrix input.")
    }

    # make counts summary
    case_counts <- unfilt_counts[, group=="case"]
    con_counts <- data.frame(unfilt_counts[, group=="control"])
    counts_summary <- data.frame(ec_names=names(case_counts),
                                 case_CPM=case_cpms,
                                 case_reads=case_counts,
                                 controls_total_reads=apply(con_counts, 1, sum))

    # write counts summary table to file for reference
    counts_out <- data.frame(counts_summary,
                             passes_CPM_threshold=keep)
    counts_out <- inner_join(counts_out, tx_ec, by = "ec_names")
    counts_out <- counts_out[, c("ec_names",
                                 "contigs",
                                 "case_reads",
                                 "controls_total_reads",
                                 "case_CPM",
                                 "passes_CPM_threshold")]
    write.table(counts_out, file=paste(outdir, "counts_summary.txt", sep="/"), row.names=FALSE, quote=FALSE, sep="\t")

    # perform DE
    n_controls <- ncol(ec_matrix) - 4
    if (n_controls == 1) {
        # add manual dispersion parameter, perform exact test
        dge <- DGEList(counts = counts, group = group)
        dge$samples <- norm_factors
        et <- exactTest(dge, dispersion=0.1)
        dx_df <- data.frame(topTags(et, n=Inf))
    } else {
        des <- model.matrix(~group)
        colnames(des)[2] <- "case"
        dge <- DGEList(counts = counts, group = group)
        dge$samples <- norm_factors
        dge <- estimateDisp(dge)
        fit <- glmQLFit(dge, design=des, robust=TRUE)
        lrt <- glmLRT(fit, coef=2)
        dx_df <- data.frame(topTags(lrt, n=Inf))
    }

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
    dx_df <- dx_df[order(dx_df$PValue, decreasing=FALSE),]

    # write full results
    significant <- dx_df$FDR < qval
    passes_logFC_threshold <- dx_df$logFC > min_logfc
    write.table(data.frame(dx_df, significant, passes_logFC_threshold),
                file=paste(outdir, "full_edgeR_results.txt", sep="/"), row.names=FALSE, quote=FALSE, sep="\t")

    dx_df <- dx_df[significant & passes_logFC_threshold,]
    return(dx_df)
}
