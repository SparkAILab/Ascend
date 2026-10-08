## ======================================================================
## realdata/drosophila/dgrp_data.R
## Builds the DGRP input for ASCEND exactly as dgrp_ascend.R does (same
## filters, same seed, same 500 SNPs and 250 genes), and records how many
## SNPs and genes survive each filter.
##
## Reviewer comments: 7 (DGRP controls: shared input for every control run),
##   minor 7 (the SNP and gene counts of the filtering flow diagram)
## How to run: sourced by dgrp_controls.R; nothing to run by hand.
##   Data (not in git, too large): put these two files in DGRP_DATA
##   (default: realdata/drosophila/):
##     dgrp2.tgeno.txt   from data.tar.gz at https://zenodo.org/records/14871341
##     GSE117850_DGRP_GEO_Table_4_Gene_Male_Line_Means.txt.gz   from GEO GSE117850
## Full write-up: docs/REVISION_REPORT.md
## ======================================================================

library(data.table)

## Read both files, filter as in the paper, cache the result as an .rds so
## the hundreds of control runs do not re-read the 1.4 GB genotype file.
dgrp_input <- function(data_dir, cache = file.path(data_dir, "dgrp_input.rds")) {
  if (file.exists(cache)) return(readRDS(cache))
  flow <- list()

  cat("Loading data...\n")
  geno <- fread(file.path(data_dir, "dgrp2.tgeno.txt"))
  expr <- fread(file.path(data_dir, "GSE117850_DGRP_GEO_Table_4_Gene_Male_Line_Means.txt.gz"))
  geno_cols <- grep("^line_", names(geno), value = TRUE)
  flow$snps_all  <- nrow(geno)
  flow$lines_geno <- length(geno_cols)

  # SNPs with no missing call ("-") in any line
  complete <- geno[, .(complete = all(!grepl("-", .SD))), .SDcols = geno_cols, by = 1:nrow(geno)]$complete
  geno <- geno[which(complete)]
  flow$snps_complete <- nrow(geno)
  for (col in geno_cols) set(geno, j = col, value = as.numeric(geno[[col]]))

  # lines with both genotype and expression
  geno_ids <- as.numeric(gsub("line_", "", geno_cols))
  expr_ids <- as.numeric(gsub("_M_mean", "", names(expr)[-1]))
  common_lines <- intersect(geno_ids, expr_ids)
  flow$lines_expr   <- length(expr_ids)
  flow$lines_common <- length(common_lines)

  # common variants, then 500 drawn at random (seed 123, as in dgrp_ascend.R)
  geno_mat <- as.matrix(geno[, paste0("line_", common_lines), with = FALSE])
  rownames(geno_mat) <- geno$snp_id
  rm(geno); gc()
  alt_freq <- rowMeans(geno_mat) / 2
  maf_vals <- pmin(alt_freq, 1 - alt_freq)
  good_maf <- maf_vals > 0.2 & maf_vals < 0.5
  flow$snps_maf <- sum(good_maf)
  if (sum(good_maf) >= 500) {
    set.seed(123)
    selected_snps <- sample(which(good_maf), 500)
  } else selected_snps <- which(good_maf)
  z_mat <- geno_mat[selected_snps, ]
  flow$snps_used <- nrow(z_mat)
  rm(geno_mat); gc()

  # genes: complete expression, mean log2 FPKM in (4, 10), variance > 0.5,
  # then the 250 with the highest variance
  expr_mat <- as.matrix(expr[, paste0(common_lines, "_M_mean"), with = FALSE])
  rownames(expr_mat) <- expr$`Gene ID`
  flow$genes_all <- nrow(expr_mat)
  complete_genes <- rowSums(is.na(expr_mat)) == 0
  flow$genes_complete <- sum(complete_genes)
  expr_complete <- expr_mat[complete_genes, ]
  gene_vars  <- apply(expr_complete, 1, var)
  gene_means <- rowMeans(expr_complete)
  flow$genes_mean_4_10 <- sum(gene_means > 4 & gene_means < 10)
  good_expr <- gene_means > 4 & gene_means < 10 & gene_vars > 0.5
  flow$genes_filtered <- sum(good_expr)
  if (sum(good_expr) >= 250) {
    top_genes <- names(sort(gene_vars[good_expr], decreasing = TRUE))[1:250]
  } else top_genes <- names(gene_vars[good_expr])
  flow$genes_used <- length(top_genes)

  out <- list(z = t(z_mat),                                   # lines x SNPs
              x = t(expr_complete[top_genes, ]),              # lines x genes
              x_pool = t(expr_complete[good_expr, ]),         # every gene passing the filter
              flow = unlist(flow))
  tmp <- paste0(cache, ".tmp", Sys.getpid())          # atomic: parallel runs never see half a file
  saveRDS(out, tmp); invisible(file.rename(tmp, cache))
  out
}

## ASCEND input from a lines x SNPs and a lines x genes matrix, standardised
## and named z1.., x1.. as in dgrp_ascend.R
dgrp_ascend_obj <- function(z, x) {
  z_std <- scale(z); x_std <- scale(x)
  dat <- cbind(z_std, x_std)
  colnames(dat) <- c(paste0("z", seq_len(ncol(z))), paste0("x", seq_len(ncol(x))))
  list(dat = as.data.frame(dat),
       metadata = list(n_samples = nrow(z), n_snps = ncol(z), n_genes = ncol(x),
                       snp_ids = colnames(z), gene_ids = colnames(x)))
}

## Edge counts and degrees from ASCEND's matrix, with the definitions of
## dgrp_ascend.R (directed = entries equal to 1, undirected = entries equal
## to 0.5, density over d(d-1)/2 pairs)
dgrp_summary <- function(M, gene_ids) {
  rownames(M) <- colnames(M) <- gene_ids
  ones <- sum(M == 1, na.rm = TRUE); halfs <- sum(M == 0.5, na.rm = TRUE)
  d <- nrow(M)
  e <- which(M == 1, arr.ind = TRUE)
  list(n_directed = ones, n_undirected = halfs,
       density = (ones + halfs) / (d * (d - 1) / 2),
       out_degree = rowSums(M == 1, na.rm = TRUE),
       in_degree  = colSums(M == 1, na.rm = TRUE),
       edges = paste(gene_ids[e[, 1]], gene_ids[e[, 2]], sep = ">"),
       stats = attr(M, "stats"))
}
