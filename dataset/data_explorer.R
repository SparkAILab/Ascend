# plots_for_paper.R
suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(factoextra)   # for PCA plotting convenience (optional)
})

dir.create("plots", showWarnings = FALSE)

# -------- Load data --------
# Option A: load ASCEND object (has combined scaled data)
ascend_obj <- readRDS("ASCEND_OBJECT.rds")   # or ASCEND_OBJECT_condensed.rds
ascend_dat_scaled <- as.data.frame(ascend_obj$dat)   # samples x variables, scaled
z_cols <- ascend_obj$z_cols
x_cols <- ascend_obj$x_cols

# Option B: load expression and methylation processed files (if you want original matrices)
expr_z <- read.csv("EXPRESSION_log2CPM_zscored.csv", stringsAsFactors = FALSE, check.names = FALSE)
meth_z <- read.csv("PROMOTER_METHYLATION_logit_zscored.csv", stringsAsFactors = FALSE, check.names = FALSE)
# convert to matrices with rownames
rownames(expr_z) <- expr_z$ensembl_gene_id; expr_z$ensembl_gene_id <- NULL
rownames(meth_z) <- meth_z$ensembl_gene_id; meth_z$ensembl_gene_id <- NULL

# If ascend_dat_scaled exists, split into methylation (z...) and expr (x...)
if(!is.null(ascend_dat_scaled)) {
  # ensure ascend_dat_scaled is data.frame with rownames = sample names
  ascend_df <- ascend_dat_scaled
  # get methylation and expr columns by prefixes
  mz <- ascend_df %>% dplyr::select(starts_with("z"))
  ex <- ascend_df %>% dplyr::select(starts_with("x"))
  # convert to matrices with genes as columns -> transpose if needed
  # here rows are samples, columns variables (genes), which we want for PCA across samples
}

# Helper colours
heat_cols <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

# -------- 1) PCA plots (expression & methylation) --------
# Expression PCA (use expr_z: rows=genes, cols=samples) -> transpose for prcomp
# -------- 1) PCA plots (expression & methylation) --------

# Helper function to prepare matrix for PCA (removes NA rows/genes)
prepare_pca_mat <- function(mat) {
  # Remove rows (genes) where *any* sample is NA.
  # We transpose later, so we remove genes that are incomplete across samples.
  mat_complete <- mat[complete.cases(mat), ] 
  if(nrow(mat_complete) == 0) stop("No complete genes left for PCA.")
  return(mat_complete)
}

# Expression PCA
if(exists("expr_z")) {
  # Expression data is likely complete since NAs were replaced with 0s and filtered earlier.
  # But we apply the filter just in case, ensuring the matrix is numeric.
  expr_mat_pca <- as.matrix(prepare_pca_mat(expr_z)) 
  expr_pca <- prcomp(t(expr_mat_pca), center = TRUE, scale. = FALSE) 
  pca_df <- data.frame(sample = rownames(expr_pca$x), PC1 = expr_pca$x[,1], PC2 = expr_pca$x[,2])
  # ... (rest of expression PCA plot code is correct)
  # ...
  ggsave("plots/PCA_expression.png", width = 7, height = 5, dpi = 200)
}

# Methylation PCA (CRITICAL FIX)
if(exists("meth_z")) {
  meth_mat_pca <- as.matrix(prepare_pca_mat(meth_z)) # Filtered matrix
  meth_pca <- prcomp(t(meth_mat_pca), center = TRUE, scale. = FALSE) # Should now work!
  pca_df <- data.frame(sample = rownames(meth_pca$x), PC1 = meth_pca$x[,1], PC2 = meth_pca$x[,2])
  # Extract variance explained for axes labels
  var_exp <- summary(meth_pca)$importance[2,]
  ggplot(pca_df, aes(PC1, PC2, label = sample)) +
    geom_point(size = 2) + geom_text_repel(size = 3, max.overlaps = 20) +
    ggtitle("PCA of promoter methylation (logit z-scored, complete genes)") + theme_minimal() +
    xlab(paste0("PC1 (", round(100 * var_exp[1],1), "%)")) +
    ylab(paste0("PC2 (", round(100 * var_exp[2],1), "%)"))
  ggsave("plots/PCA_methylation.png", width = 7, height = 5, dpi = 200)
}

# -------- 2) Heatmap: top variable genes (expression) --------
topN <- 100
if(exists("expr_z")) {
  gene_vars <- apply(expr_z, 1, var, na.rm = TRUE)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(topN, length(gene_vars))]
  mat <- expr_z[top_genes, , drop = FALSE]  # genes x samples
  
  # FIX: Remove any gene (row) that has NA in the topN subset
  mat_complete <- mat[complete.cases(mat), ] 
  
  if(nrow(mat_complete) > 1) { # Need at least 2 rows for clustering
    pheatmap(mat_complete, color = heat_cols, show_rownames = FALSE, show_colnames = FALSE,
             main = paste0("Expression heatmap: top ", nrow(mat_complete), " complete var genes"),
             filename = "plots/heatmap_expression_topvars.pdf", width = 8, height = 10)
  } else {
    message("Not enough complete genes for Expression Heatmap clustering.")
  }
}

# -------- 3) Heatmap: top variable promoters (methylation) --------
# -------- 3) Heatmap: top variable promoters (methylation) --------
if(exists("meth_z")) {
  gene_vars_m <- apply(meth_z, 1, var, na.rm = TRUE)
  top_meth_genes <- names(sort(gene_vars_m, decreasing = TRUE))[1:min(topN, length(gene_vars_m))]
  mat_m <- meth_z[top_meth_genes, , drop = FALSE]
  
  # CRITICAL FIX: Use only complete rows (genes) for clustering
  mat_m_complete <- mat_m[complete.cases(mat_m), ] 
  
  if(nrow(mat_m_complete) > 1) { # Need at least 2 rows for clustering
    pheatmap(mat_m_complete, color = heat_cols, show_rownames = FALSE, show_colnames = FALSE,
             main = paste0("Methylation heatmap: top ", nrow(mat_m_complete), " complete var promoters"),
             filename = "plots/heatmap_methylation_topvars.pdf", width = 8, height = 10)
  } else {
    message("Not enough complete genes for Methylation Heatmap clustering.")
  }
}

# -------- 4) Distribution of methylation beta values --------
if(file.exists("PROMOTER_METHYLATION_beta.csv")) {
  bet <- read.csv("PROMOTER_METHYLATION_beta.csv", stringsAsFactors = FALSE, check.names = FALSE)
  bet_long <- bet %>% pivot_longer(-ensembl_gene_id, names_to = "sample", values_to = "beta")
  ggplot(bet_long, aes(x = beta)) + geom_histogram(bins = 80) +
    facet_wrap(vars(sample), scales = "free_y", ncol = 5) +
    ggtitle("Beta-value distributions per sample (promoter methylation)") +
    xlab("Beta (0..1)") + theme_minimal()
  ggsave("plots/methylation_beta_hist_by_sample.png", width = 12, height = 10, dpi = 200)
  # overall density
  ggplot(bet_long, aes(x = beta)) + geom_density(na.rm = TRUE) +
    ggtitle("Overall distribution of promoter beta-values") + xlab("Beta (0..1)") + theme_minimal()
  ggsave("plots/methylation_beta_density.png", width = 6, height = 4, dpi = 200)
}

# -------- 5) Paired scatter plots (methylation vs expression) for select genes --------
# pick a few genes by symbol or Ensembl id. Replace with your 6-7 gene IDs.
genes_of_interest <- c("ENSRNOG00000000009", "ENSRNOG00000000012")[1:2]  # replace with your list
# create a combined long table for plotting
if(exists("expr_z") && exists("meth_z")) {
  for(g in genes_of_interest) {
    if(!(g %in% rownames(expr_z)) || !(g %in% rownames(meth_z))) {
      message("Gene not found in both matrices: ", g); next
    }
    expr_vals <- expr_z[g, , drop = FALSE] %>% t() %>% as.data.frame()
    meth_vals <- meth_z[g, , drop = FALSE] %>% t() %>% as.data.frame()
    colnames(expr_vals) <- "expr_z"
    colnames(meth_vals) <- "meth_z"
    df <- bind_cols(expr_vals, meth_vals) %>% mutate(sample = rownames(expr_vals))
    # scatter with linear fit and Pearson r
    cor_val <- cor(df$expr_z, df$meth_z, use = "complete.obs")
    p <- ggplot(df, aes(x = meth_z, y = expr_z)) +
      geom_point() + geom_smooth(method = "lm", se = TRUE) +
      ggtitle(paste0("Gene ", g, " : methylation vs expression (z-scored)\nPearson r = ", round(cor_val, 3))) +
      theme_minimal()
    ggsave(filename = paste0("plots/paired_scatter_", g, ".png"), plot = p, width = 6, height = 5, dpi = 200)
  }
}

# -------- 6) Sample-sample correlation heatmap (expression) --------
if(exists("expr_z")) {
  samp_cor <- cor(expr_z, use = "pairwise.complete.obs")  # cor between samples (columns)
  pheatmap(samp_cor, main = "Sample-sample correlation (expression)", filename = "plots/sample_correlation_expression.pdf",
           color = colorRampPalette(brewer.pal(9,"Blues"))(100), width = 8, height = 8)
}


# Note: This code assumes expr_z and meth_z are defined and are gene (row) x sample (column) matrices.

# -------- 7) Global Scatter Plot: Methylation vs. Expression (QC) --------
# Purpose: Visually assess the overall correlation between the two Z-scored datasets.
message("Generating Global Scatter Plot: Methylation vs. Expression...")

if (exists("expr_z") && exists("meth_z")) {
  # 1. Convert matrices to long format
  ex_long <- expr_z %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "expr_z")
  
  mz_long <- meth_z %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "meth_z")
  
  # 2. Join the data on gene and sample to ensure paired observations
  combined_long <- inner_join(ex_long, mz_long, by = c("gene", "sample"))
  
  # 3. Compute correlation over all paired observations
  global_cor <- cor(combined_long$expr_z, combined_long$meth_z, use = "complete.obs")
  
  # 4. Plot (sampling 10% for efficiency and clarity due to high density)
  p <- ggplot(combined_long %>% sample_frac(0.1), 
              aes(x = meth_z, y = expr_z)) +
    geom_point(alpha = 0.2, color = "darkgrey") +
    # 2D density contours show where the majority of the data lies
    geom_density_2d(colour = "#1f78b4") + 
    # Linear fit shows the trend
    geom_smooth(method = "lm", se = FALSE, color = "#e31a1c", linewidth = 1.2) +
    ggtitle(paste0("Global Methylation vs. Expression (Z-scored)\nAll Paired Gene-Samples (Pearson r = ", round(global_cor, 3), ")")) +
    xlab("Promoter Methylation (Logit Z-score)") + 
    ylab("Gene Expression (log2CPM Z-score)") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10))
  
  ggsave("plots/global_scatter_meth_vs_expr.png", plot = p, width = 7, height = 6, dpi = 200)
  message("Wrote plots/global_scatter_meth_vs_expr.png")
}


# -------- 8) Gene Density vs. Variance Plot (QC for filtering bias) --------
# Purpose: Check if variance filtering (Section 12) disproportionately removes low-expressed genes.
message("Generating Mean-Variance Plot...")

if (exists("expr_z") && exists("meth_z")) {
  # Calculate Mean and Variance for Expression
  expr_mean <- rowMeans(expr_z, na.rm = TRUE)
  expr_var <- apply(expr_z, 1, var, na.rm = TRUE)
  
  # Calculate Mean and Variance for Methylation
  meth_mean <- rowMeans(meth_z, na.rm = TRUE)
  meth_var <- apply(meth_z, 1, var, na.rm = TRUE)
  
  # Gene Filtering Status (from Section 12 logic)
  # We must use the intersection of genes available in both expr_z and meth_z
  common_genes <- intersect(names(expr_var), names(meth_var))
  
  # Create the data frame for the plot
  mean_var_df <- data.frame(
    gene = common_genes,
    expr_mean = expr_mean[common_genes],
    expr_var = expr_var[common_genes],
    meth_mean = meth_mean[common_genes],
    meth_var = meth_var[common_genes]
  ) %>%
    # Add filtering status (assuming your thresholds are var_threshold_expr = 0.2 and var_threshold_meth = 0.02)
    mutate(
      is_kept_expr = expr_var >= var_threshold_expr,
      is_kept_meth = meth_var >= var_threshold_meth,
      is_kept_both = is_kept_expr & is_kept_meth,
      status = case_when(
        is_kept_both ~ "Kept (Passed Both Filters)",
        is_kept_expr ~ "Removed (Failed Meth Var)",
        is_kept_meth ~ "Removed (Failed Expr Var)",
        TRUE ~ "Removed (Failed Both Var)"
      )
    )
  
  # Plot 8a: Expression Mean vs. Variance
  p_expr <- ggplot(mean_var_df, aes(x = expr_mean, y = expr_var, color = status)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_y_log10() + # Variance is typically log-scaled
    geom_hline(yintercept = var_threshold_expr, linetype = "dashed", color = "black") +
    labs(title = "Expression: Mean vs. Variance (Z-scored)",
         x = "Gene Mean Expression (Z-score)",
         y = "Gene Variance (log10)",
         color = "Filtering Status") +
    theme_minimal()
  ggsave("plots/mean_variance_expression.png", plot = p_expr, width = 8, height = 6, dpi = 200)
  
  # Plot 8b: Methylation Mean vs. Variance
  p_meth <- ggplot(mean_var_df, aes(x = meth_mean, y = meth_var, color = status)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_y_log10() +
    geom_hline(yintercept = var_threshold_meth, linetype = "dashed", color = "black") +
    labs(title = "Methylation: Mean vs. Variance (Z-scored)",
         x = "Promoter Mean Methylation (Z-score)",
         y = "Promoter Variance (log10)",
         color = "Filtering Status") +
    theme_minimal()
  ggsave("plots/mean_variance_methylation.png", plot = p_meth, width = 8, height = 6, dpi = 200)
  
  message("Wrote plots/mean_variance_expression.png and plots/mean_variance_methylation.png")
}













# -------- Done --------
message("Plots written to ./plots/ — check files and adjust gene list / topN as needed.")
