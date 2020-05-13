library(dplyr)
library(data.table)

get_significant_ecs <- function(case_name, ec_matrix, tx_ec, outdir, minCaseCount=minCaseCount, n_zero_controls=n_zero_controls) {
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
    keep <- as.numeric(counts[, group=="case"]) >= minCaseCount
    counts <- counts[keep,]
    if(nrow(counts) == 0) {
        stop(paste("All ECCs were below the min case count cutoff of", paste0(minCaseCount, ".")))
    }
    counts <- counts[rownames(counts)%in%tx_ec$ec_names,]
    if(nrow(data) == 0) {
        stop("No novel ECs exist in EC matrix. Please check EC matrix input.")
    }

    # make counts summary
    case_counts <- counts[, group=="case"]
    con_counts <- data.frame(counts[, group=="control"])
    counts_summary <- data.frame(ec_names=names(case_counts),
                                 case_reads=case_counts,
                                 controls_total_reads=apply(con_counts, 1, sum))

    # write counts table to file for reference
    counts_out <- data.frame(ec_names=rownames(counts), counts)
    write.table(counts_out, file=paste(outdir, "counts_table.txt", sep="/"), row.names=FALSE, quote=FALSE, sep="\t")

    # perform DE
    if (test) {
        # add manual dispersion parameter, perform exact test
        dge <- DGEList(counts = counts, group = group)
        dge$samples <- norm_factors
        et <- exactTest(dge, dispersion = 0.1)
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
    print(paste("There were", sum(significant), "significant novel-contig ECs."))

    # get number of contigs in each EC
    contigs_in_ec <- data.frame(table(tx_ec$ec_names))
    contigs_in_ec$Var1 <- as.character(contigs_in_ec$Var1)
    colnames(contigs_in_ec) <- c("ec_names", "n_contigs_in_ec")

    dx_df <- contigs_in_ec[contigs_in_ec$ec_names%in%sig_ecs,]
    dx_df <- left_join(dx_df, tx_ec, by="ec_names")
    dx_df <- left_join(dx_df, counts_summary, by="ec_names")

    return(dx_df)
}
