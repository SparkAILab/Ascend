# Assuming expr_final and meth_final are loaded and aligned (gene x sample)
expr_z <- read.csv("EXPRESSION_log2CPM_zscored.csv", stringsAsFactors = FALSE, check.names = FALSE)
meth_z <- read.csv("PROMOTER_METHYLATION_logit_zscored.csv", stringsAsFactors = FALSE, check.names = FALSE)
# convert to matrices with rownames


# Only consider genes with sufficient methylation data
meth_nonNA_counts <- rowSums(!is.na(meth_z))
keep_presence <- meth_nonNA_counts >= min_samples_with_meth

# Calculate variance (handling NAs)
expr_sd <- apply(expr_z, 1, sd, na.rm = TRUE)
meth_sd <- apply(meth_z, 1, sd, na.rm = TRUE)
keep_variance <- (expr_sd >= var_threshold_expr) & (meth_sd >= var_threshold_meth)

keep_genes <- keep_presence & keep_variance
kept_gene_ids <- rownames(expr_z)[keep_genes]
message("Genes kept after filtering: ", length(kept_gene_ids))





expr_final <- expr_z[kept_gene_ids, , drop = FALSE]
meth_final <- meth_z[kept_gene_ids, , drop = FALSE]
# 1. Align and Filter to ensure only common genes are analyzed
common_genes <- intersect(rownames(expr_final), rownames(meth_final))
expr_sub <- expr_final[common_genes, ]
meth_sub <- meth_final[common_genes, ]

# 2. Calculate Pearson correlation for each gene across all samples
gene_correlations <- apply(expr_sub, 1, function(expr_row) {
  gene_id <- rownames(expr_sub)[which(apply(expr_sub, 1, identical, expr_row))]
  meth_row <- meth_sub[gene_id, ]
  
  # Ensure the row exists and is aligned
  if(length(meth_row) == 0) return(NA)
  
  # Use pairwise complete observation to handle NAs in methylation
  cor_val <- cor(as.numeric(expr_row), as.numeric(meth_row), use = "pairwise.complete.obs")
  return(cor_val)
})

# 3. Convert to data frame and sort
cor_df <- data.frame(
  ensembl_gene_id = names(gene_correlations),
  r_value = gene_correlations
) %>% 
  filter(!is.na(r_value)) %>%
  arrange(r_value) # Sort in ascending order (most negative first)

# 4. View Top 10 Anti-correlated Candidates for X Variables
print(head(cor_df, 10))

# Now manually cross-reference these 10 genes with your biological tables (like image_147e48.png)