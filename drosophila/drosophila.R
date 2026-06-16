# data sources: Expression: GSE117850 on GEO
# Genotype: https://zenodo.org/records/14871341 \\ scroll down and download the data.tar.gz which includes the necessary geno files.





# ============================================
# STEP 1: MEMORY-EFFICIENT DATA LOADING
# ============================================
library(data.table)
library(Matrix)

# 1. LOAD AND CONVERT (keeping as character for now)
cat("Loading data...\n")
geno <- fread("dgrp2.tgeno.txt")
expr <- fread("GSE117850_DGRP_GEO_Table_4_Gene_Male_Line_Means.txt.gz")

# 2. IDENTIFY COMPLETE CASES
# Find SNPs with NO missing data across all lines
geno_cols <- grep("^line_", names(geno), value = TRUE)

# Fast: Count NAs per SNP (still as character)
complete_snps <- geno[, .(complete = all(!grepl("-", .SD))), 
                      .SDcols = geno_cols, by = 1:nrow(geno)]

complete_snp_indices <- which(complete_snps$complete)
cat("\nSNPs with NO missing data:", length(complete_snp_indices), "\n")

# 3. FILTER TO COMPLETE SNPS AND CONVERT
geno_complete <- geno[complete_snp_indices]
rm(geno); gc()

# Convert to numeric (0/2) - ALL values are now valid
for (col in geno_cols) {
  set(geno_complete, j = col, value = as.numeric(geno_complete[[col]]))
}

# 4. MATCH SAMPLES
geno_ids <- as.numeric(gsub("line_", "", geno_cols))
expr_ids <- as.numeric(gsub("_M_mean", "", names(expr)[-1]))

common_lines <- intersect(geno_ids, expr_ids)
cat("Common lines:", length(common_lines), "\n")

# 5. SELECT OPTIMAL 500 SNPS FROM COMPLETE DATA
geno_mat <- as.matrix(geno_complete[, paste0("line_", common_lines), with = FALSE])
rownames(geno_mat) <- geno_complete$snp_id

# Calculate MAF (no NAs!)
alt_freq <- rowMeans(geno_mat) / 2  # All 0 or 2
maf_vals <- pmin(alt_freq, 1 - alt_freq)

# Select SNPs with good MAF
good_maf <- maf_vals > 0.2 & maf_vals < 0.5
cat("SNPs with MAF 0.2-0.5:", sum(good_maf), "\n")

# If we have enough, sample 500
if(sum(good_maf) >= 500) {
  set.seed(123)
  selected_snps <- sample(which(good_maf), 500)
} else {
  # Take all good MAF SNPs
  selected_snps <- which(good_maf)
  cat("Warning: Only", length(selected_snps), "SNPs meet criteria\n")
}

z_mat <- geno_mat[selected_snps, ]

# 6. SELECT GENES WITH COMPLETE EXPRESSION
expr_mat <- as.matrix(expr[, paste0(common_lines, "_M_mean"), with = FALSE])
rownames(expr_mat) <- expr$`Gene ID`

# Find genes with NO missing expression
complete_genes <- rowSums(is.na(expr_mat)) == 0
cat("Genes with NO missing expression:", sum(complete_genes), "\n")

# Filter and select top variable genes
expr_complete <- expr_mat[complete_genes, ]
gene_vars <- apply(expr_complete, 1, var)
gene_means <- rowMeans(expr_complete)

# Select genes with good expression levels
good_expr <- gene_means > 4 & gene_means < 10 & gene_vars > 0.5
cat("Genes with good expression levels:", sum(good_expr), "\n")

# Sort by variance and take top 250
if(sum(good_expr) >= 250) {
  top_genes <- names(sort(gene_vars[good_expr], decreasing = TRUE))[1:250]
} else {
  top_genes <- names(gene_vars[good_expr])
  cat("Warning: Only", length(top_genes), "genes meet criteria\n")
}

x_mat <- expr_complete[top_genes, ]

# 7. FINAL NA-FREE DATASET
cat("\n=== FINAL DATASET (NO NAS) ===\n")
cat("Samples:", ncol(z_mat), "\n")
cat("Z (SNPs):", nrow(z_mat), "\n") 
cat("X (Genes):", nrow(x_mat), "\n")
cat("Missing values in Z:", sum(is.na(z_mat)), "\n")
cat("Missing values in X:", sum(is.na(x_mat)), "\n")

# Standardize (no NAs!)
z_std <- scale(t(z_mat))  # samples × SNPs
x_std <- scale(t(x_mat))  # samples × genes

# 8. CREATE ASCEND INPUT
ascend_data <- cbind(z_std, x_std)
colnames(ascend_data) <- c(
  paste0("z", 1:ncol(z_std)),
  paste0("x", 1:ncol(x_std))
)

real_obj <- list(
  dat = as.data.frame(ascend_data),
  metadata = list(
    n_samples = ncol(z_mat),
    n_snps = nrow(z_mat),
    n_genes = nrow(x_mat),
    snp_ids = rownames(z_mat),
    gene_ids = rownames(x_mat)
  )
)

# 9. RUN ASCEND (CLEAN DATA)
source("ascend_fn.R")

cat("\nRunning ASCEND on CLEAN data...\n")
result <- ascend_fn(real_obj, alpha = 0.05, fdr_correction = TRUE, maxiter = 10)




#results saved as real_obj_result.rds
#=============================================

# ============================================
# FIXED ASCEND RESULTS ANALYSIS
# ============================================
library(tidyverse)

# Load your results
result <- real_obj_result  # Your ASCEND adjacency matrix

# ============================================
# 1. BASIC STATISTICS (FIXED)
# ============================================
cat("=== ASCEND RESULTS SUMMARY ===\n")
cat("Matrix dimensions:", dim(result), "\n\n")

# Check if result has row/column names
if(is.null(rownames(result))) {
  cat("Warning: No row names in result matrix\n")
  rownames(result) <- paste0("gene_", 1:nrow(result))
}
if(is.null(colnames(result))) {
  cat("Warning: No column names in result matrix\n")
  colnames(result) <- paste0("gene_", 1:ncol(result))
}

# Count different edge types
zeros <- sum(result == 0, na.rm = TRUE)
ones <- sum(result == 1, na.rm = TRUE)
halfs <- sum(result == 0.5, na.rm = TRUE)
nas <- sum(is.na(result))

cat("Edge counts:\n")
cat("  0 (no edge):", zeros, "\n")
cat("  1 (directed):", ones, "\n")
cat("  0.5 (undirected):", halfs, "\n")
cat("  NA (untested):", nas, "\n\n")

# Network density
total_edges <- ones + halfs
possible_edges <- nrow(result) * (ncol(result) - 1) / 2
density <- total_edges / possible_edges
cat("Network density:", round(density * 100, 2), "%\n")

# ============================================
# 2. EXTRACT EDGES (ROBUST VERSION)
# ============================================

# Function to safely extract edges
extract_edges <- function(adj_matrix, value) {
  edges <- which(adj_matrix == value, arr.ind = TRUE)
  
  if(length(edges) == 0) {
    return(data.frame(from = character(0), to = character(0)))
  }
  
  # Ensure edges is a matrix (not a vector)
  if(!is.matrix(edges)) {
    edges <- matrix(edges, nrow = 1)
  }
  
  data.frame(
    from = rownames(adj_matrix)[edges[, 1]],
    to = colnames(adj_matrix)[edges[, 2]],
    value = value
  )
}

# Extract all edge types
directed_edges <- extract_edges(result, 1)
undirected_edges <- extract_edges(result, 0.5)

cat("\n=== DIRECTED EDGES ===\n")
if(nrow(directed_edges) > 0) {
  cat("Found", nrow(directed_edges), "directed edges\n")
  print(head(directed_edges, 10))
  
  # Save to file
  write.csv(directed_edges, "ascend_directed_edges.csv", row.names = FALSE)
  cat("\nSaved to: ascend_directed_edges.csv\n")
} else {
  cat("No directed edges found\n")
}

cat("\n=== UNDIRECTED EDGES ===\n")
if(nrow(undirected_edges) > 0) {
  cat("Found", nrow(undirected_edges), "undirected edges\n")
  print(head(undirected_edges, 10))
  
  write.csv(undirected_edges, "ascend_undirected_edges.csv", row.names = FALSE)
  cat("\nSaved to: ascend_undirected_edges.csv\n")
} else {
  cat("No undirected edges found\n")
}

# ============================================
# 3. NODE STATISTICS
# ============================================
cat("\n=== NODE STATISTICS ===\n")

# Calculate indegree and outdegree safely
calculate_degree <- function(adj_matrix) {
  gene_names <- rownames(adj_matrix)
  
  indegree <- setNames(rep(0, length(gene_names)), gene_names)
  outdegree <- setNames(rep(0, length(gene_names)), gene_names)
  
  if(nrow(directed_edges) > 0) {
    # Count outgoing edges (from → to)
    out_counts <- table(directed_edges$from)
    outdegree[names(out_counts)] <- out_counts
    
    # Count incoming edges
    in_counts <- table(directed_edges$to)
    indegree[names(in_counts)] <- in_counts
  }
  
  data.frame(
    gene = gene_names,
    indegree = as.numeric(indegree),
    outdegree = as.numeric(outdegree),
    total_edges = as.numeric(indegree + outdegree)
  )
}

node_stats <- calculate_degree(result)

# Show top nodes
cat("\nTop 10 genes by total edges:\n")
print(head(node_stats[order(-node_stats$total_edges), ], 10))

# Save node statistics
write.csv(node_stats, "ascend_node_statistics.csv", row.names = FALSE)
cat("\nNode statistics saved: ascend_node_statistics.csv\n")

# ============================================
# 4. QUICK VISUALIZATION (SIMPLE)
# ============================================
if(nrow(directed_edges) > 0) {
  cat("\n=== NETWORK VISUALIZATION ===\n")
  
  # Simple text-based visualization
  cat("\nNetwork structure (first 20 edges):\n")
  for(i in 1:min(20, nrow(directed_edges))) {
    cat(directed_edges$from[i], "→", directed_edges$to[i], "\n")
  }
  
  # Simple plot if we have edges
  if(nrow(directed_edges) <= 50) {  # Only plot small networks
    tryCatch({
      library(igraph)
      
      # Create graph
      g <- graph_from_data_frame(directed_edges[, 1:2], directed = TRUE)
      
      # Basic plot
      pdf("ascend_simple_network.pdf", width = 8, height = 6)
      plot(g, 
           vertex.size = 10,
           vertex.color = "lightblue",
           vertex.label.cex = 0.8,
           edge.arrow.size = 0.5,
           layout = layout_with_fr,
           main = paste("ASCEND Network (", nrow(directed_edges), "edges)"))
      dev.off()
      cat("\nNetwork plot saved: ascend_simple_network.pdf\n")
    }, error = function(e) {
      cat("Could not create network plot:", e$message, "\n")
    })
  }
}

# ============================================
# 5. MATRIX INSPECTION
# ============================================
cat("\n=== MATRIX INSPECTION ===\n")

# Show a subset of the matrix
cat("\nFirst 10x10 submatrix:\n")
print(result[1:min(10, nrow(result)), 1:min(10, ncol(result))])

# Check for patterns
cat("\nMatrix pattern:\n")
cat("- Upper triangle has", sum(result[upper.tri(result)] == 1, na.rm = TRUE), "ones\n")
cat("- Lower triangle has", sum(result[lower.tri(result)] == 1, na.rm = TRUE), "ones\n")

# ============================================
# 6. DIAGNOSTIC CHECK
# ============================================
cat("\n=== DIAGNOSTIC CHECK ===\n")

# Check if matrix is upper/lower triangular
is_upper <- all(result[lower.tri(result)] == 0, na.rm = TRUE)
is_lower <- all(result[upper.tri(result)] == 0, na.rm = TRUE)

cat("Matrix is upper triangular:", is_upper, "\n")
cat("Matrix is lower triangular:", is_lower, "\n")

# Check DAG property
tryCatch({
  binary_adj <- (result == 1 | result == 0.5) * 1
  binary_adj[is.na(binary_adj)] <- 0
  
  g_test <- graph_from_adjacency_matrix(binary_adj)
  is_dag <- igraph::is_dag(g_test)
  cat("Is DAG (no cycles):", is_dag, "\n")
}, error = function(e) {
  cat("Could not check DAG property\n")
})

# ============================================
# 7. SAVE FULL REPORT
# ============================================
sink("ascend_analysis_report.txt")
cat("ASCEND Analysis Report\n")
cat("Generated:", Sys.time(), "\n\n")
cat("Matrix dimensions:", dim(result), "\n")
cat("Directed edges:", nrow(directed_edges), "\n")
cat("Undirected edges:", nrow(undirected_edges), "\n")
cat("Network density:", round(density * 100, 2), "%\n\n")

if(nrow(directed_edges) > 0) {
  cat("Top regulators:\n")
  top_regulators <- node_stats[order(-node_stats$outdegree), ][1:5, ]
  print(top_regulators)
}
sink()

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Full report saved: ascend_analysis_report.txt\n")
