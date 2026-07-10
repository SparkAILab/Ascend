# data sources: Expression: GSE117850 on GEO and download GSE117850_DGRP_GEO_Table_4_Gene_Male_Line_Means.txt.gz
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
#source("ascend.R")

cat("\nRunning ASCEND on CLEAN data...\n")
result <- ascend(real_obj, alpha = 0.05, alpha_mb = 0.05, fdr = TRUE, min_votes = 1)




#results saved as real_obj_result.rds
#=============================================

# ============================================
# FIXED ASCEND RESULTS ANALYSIS
# ============================================
library(tidyverse)



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






# ======================================================================
#  DYNAMICALLY IDENTIFY HUBS (INDEX-CORRECTED)
# ======================================================================
cat("\n=== EXTRACTING EMPIRICAL HUBS WITH REAL FLYBASE IDs ===\n")

# 1. Compute strict causal out-degrees
causal_out_degrees <- rowSums(result == 1, na.rm = TRUE)

# 2. Extract the index number from placeholders (e.g., "x170" -> 170)
generic_names <- names(causal_out_degrees)
gene_indices  <- as.numeric(gsub("x", "", generic_names))

# 3. Map indices back to the true FlyBase IDs stored in metadata
true_flybase_ids <- real_obj$metadata$gene_ids[gene_indices]

# 4. Build the clean hub dataframe
df_hubs <- data.frame(
  generic_id = generic_names,
  gene_id    = true_flybase_ids, # This now contains actual FBgn IDs!
  out_degree = as.numeric(causal_out_degrees),
  stringsAsFactors = FALSE
)

# Sort to isolate the top drivers
df_hubs <- df_hubs[order(-df_hubs$out_degree), ]

cat("Top discovered regulators mapped to FlyBase IDs:\n")
print(head(df_hubs, 10))

# ======================================================================
#  GENOMIC CHROMOSOME MAPPING VIA BIOMART
# ======================================================================
cat("\nConnecting to Ensembl BioMart for chromosome coordinates...\n")
library(biomaRt)
library(ggplot2)
library(dplyr)

tryCatch({
  mart <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")
  
  # Query using genuine FlyBase IDs
  coords <- getBM(
    attributes = c("flybase_gene_id", "external_gene_name", "chromosome_name", "start_position"),
    filters    = "flybase_gene_id",
    values     = df_hubs$gene_id,
    mart       = mart
  )
  
  colnames(coords) <- c("gene_id", "gene_symbol", "chromosome", "position")
  df_genomic_hubs <- merge(df_hubs, coords, by = "gene_id", all.x = TRUE)
  
  # Fallback for empty symbols
  df_genomic_hubs$gene_symbol[is.na(df_genomic_hubs$gene_symbol) | df_genomic_hubs$gene_symbol == ""] <- 
    df_genomic_hubs$gene_id[is.na(df_genomic_hubs$gene_symbol) | df_genomic_hubs$gene_symbol == ""]
  
  # Filter out scaffolds, focus on major chromosomal arms
  df_genomic_hubs <- df_genomic_hubs %>%
    filter(chromosome %in% c("2L", "2R", "3L", "3R", "X", "4")) %>%
    mutate(chromosome = factor(chromosome, levels = c("X", "2L", "2R", "3L", "3R", "4"))) %>%
    arrange(chromosome, position)
  
  write.csv(df_genomic_hubs, "ascend_discovered_hubs.csv", row.names = FALSE)
  cat("Genomic positions mapped! Saved to 'ascend_discovered_hubs.csv'\n")
  
  cat("\nTop Chromosomal Hubs:\n")
  print(head(df_genomic_hubs %>% arrange(-out_degree), 10))
  
  # ======================================================================
  #  GENERATE CAUSAL MANHATTAN PLOT
  # ======================================================================
  cat("Generating Manhattan plot (manhattan.pdf)...\n")
  
  pdf("manhattan.pdf", width = 11, height = 5)
  p <- ggplot(df_genomic_hubs, aes(x = position / 1e6, y = out_degree, color = chromosome)) +
    geom_point(aes(size = out_degree), alpha = 0.8) +
    geom_segment(aes(x = position / 1e6, xend = position / 1e6, y = 0, yend = out_degree), alpha = 0.4) +
    facet_grid(. ~ chromosome, scales = "free_x", space = "free_x") +
    scale_color_brewer(palette = "Set2") +
    labs(
      title = "Genomic Architecture of Discovered Causal Influence (ASCEND)",
      x = "Genomic Position (Mb)",
      y = "Causal Out-Degree (Targets Regulated)",
      size = "Out-degree"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      panel.border = element_rect(color = "grey85", fill = NA),
      strip.background = element_rect(fill = "grey95", color = "grey85"),
      panel.grid.minor = element_blank()
    )
  print(p)
  dev.off()
  cat("Plot successfully generated and saved to 'manhattan.pdf'\n")
  
}, error = function(e) {
  cat("BioMart connection failed. Error message:", e$message, "\n")
  cat("Saving non-mapped empirical results to 'ascend_discovered_hubs_raw.csv'\n")
  write.csv(df_hubs, "ascend_discovered_hubs_raw.csv", row.names = FALSE)
})














# ======================================================================
# GENERATE CLASSIC SINGLE-PANEL GENOMIC MAP (DYNAMIC RANGE FIX)
# ======================================================================
cat("Generating classic style genomic map...\n")

library(ggplot2)
library(ggrepel)
library(dplyr)

# 1. Prepare continuous coordinates across combined chromosomes
df_classic <- df_genomic_hubs %>%
  filter(chromosome %in% c("X", "2L", "2R", "3L", "3R")) %>%
  mutate(chromosome = factor(chromosome, levels = c("X", "2L", "2R", "3L", "3R"))) %>%
  group_by(chromosome) %>%
  arrange(position) %>%
  ungroup()

# Create a continuous index for the X-axis so arms sit side-by-side
chrom_lengths <- df_classic %>% 
  group_by(chromosome) %>% 
  summarise(max_pos = max(position, na.rm = TRUE)) %>% 
  mutate(offset = lag(cumsum(as.numeric(max_pos)), default = 0))

df_classic <- df_classic %>%
  inner_join(chrom_lengths, by = "chromosome") %>%
  mutate(global_position = (position + offset) / 1e6)

# Calculate midpoints for the X-axis chromosome labels
chrom_labels <- df_classic %>%
  group_by(chromosome) %>%
  summarise(center = mean(global_position))

# 2. DYNAMICALLY DEFINE THRESHOLDS BASED ON DATA DETECTED
max_out_degree <- max(df_classic$out_degree, na.rm = TRUE)

# If it's a high-degree run, label top 10 hubs; if it's lower, label those >= 5
if (max_out_degree > 10) {
  df_labels <- df_classic %>% arrange(-out_degree) %>% head(10)
  y_max_limit <- max_out_degree * 1.1 # Give 10% breathing room at the top
} else {
  df_labels <- df_classic %>% filter(gene_symbol %in% c("CG15034", "CG17928", "CG13323", "CG6839", "Gba1a") | out_degree >= 5)
  y_max_limit <- 8.5
}

# 3. Plotting using the clean base-R aesthetic
pdf("manhattan_classic_style.pdf", width = 9, height = 7)

p_classic <- ggplot(df_classic, aes(x = global_position, y = out_degree, color = chromosome)) +
  geom_point(size = 3, alpha = 0.8) +
  # Label the specific master regulators cleanly without overlaps
  geom_text_repel(
    data = df_labels, 
    aes(label = gene_symbol),
    color = "black", 
    size = 3.5,
    box.padding = 0.6,
    point.padding = 0.4,
    max.overlaps = 15,
    fontface = "plain",
    family = "serif"
  ) +
  scale_x_continuous(breaks = chrom_labels$center, labels = chrom_labels$chromosome) +
  scale_y_continuous(limits = c(0, y_max_limit)) +
  # Custom color palette matching the distinct blocks in your target PDF
  scale_color_manual(values = c("X" = "#E41A1C", "2L" = "#377EB8", "2R" = "#4DAF4A", "3L" = "#984EA3", "3R" = "#FF7F00")) +
  labs(
    x = "Chromosome",
    y = "Causal Out-degree (Number of Targets)"
  ) +
  theme_classic() + 
  theme(
    legend.position = "none",
    text = element_text(family = "serif", size = 14),
    # Fixed deprecated 'size' to 'linewidth'
    axis.line = element_line(linewidth = 0.6, color = "black"),
    axis.ticks = element_line(linewidth = 0.6, color = "black"),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
  )

print(p_classic)
dev.off()

cat("Classic style plot successfully refreshed with all points captured!\n")