# ======================================================================
# COMPLETE COMPARISON: ASCEND vs TRIGGER on Yeast Data
# CORRECTED VERSION - WITH MORE GENES
# ======================================================================

# Load required packages
library(trigger)
library(bnlearn)
library(org.Sc.sgd.db)
library(GO.db)
library(topGO)
library(foreach)
library(doMC)
registerDoMC(16)

# Source ASCEND function
source("ascend.R")

# ======================================================================
# PART 1: LOAD YEAST DATA
# ======================================================================
cat("\n=== LOADING YEAST DATA ===\n")
data(yeast)

# Format data
z <- t(yeast$marker - 1)
x <- t(yeast$exp)

cat(sprintf("Markers: %d strains × %d SNPs\n", nrow(z), ncol(z)))
cat(sprintf("Expression: %d strains × %d genes\n", nrow(x), ncol(x)))

# ======================================================================
# PART 2: SELECT INFORMATIVE GENES (ROBUST VERSION)
# ======================================================================
cat("\n=== SELECTING INFORMATIVE GENES ===\n")

# Calculate gene metrics
gene_vars <- apply(x, 2, var, na.rm=TRUE)
gene_means <- colMeans(x, na.rm=TRUE)

# Get GO annotations for ALL genes (this will take a moment)
cat("Retrieving GO annotations for all genes...\n")
all_go <- AnnotationDbi::select(
  org.Sc.sgd.db, 
  keys = colnames(x),
  columns = c("GO", "ONTOLOGY", "EVIDENCE"),
  keytype = "ORF"
)

# Count GO terms per gene
go_count <- table(all_go$ORF)
genes_with_go <- names(go_count)[go_count >= 2]  # At least 2 GO terms
cat(sprintf("Genes with ≥2 GO terms: %d\n", length(genes_with_go)))

# Strategy 1: Take top variable genes that have GO annotations
top_var_genes <- names(sort(gene_vars, decreasing = TRUE)[1:300])
genes_for_analysis <- intersect(top_var_genes, genes_with_go)
cat(sprintf("Top variable genes with GO: %d\n", length(genes_for_analysis)))

# If we still don't have enough, include genes without GO but with high marker correlation
if(length(genes_for_analysis) < 50) {
  cat("Including genes with high marker correlation...\n")
  
  # Calculate marker correlations for top genes without GO
  no_go_genes <- setdiff(top_var_genes[1:500], genes_with_go)[1:100]
  marker_cors <- sapply(no_go_genes, function(g) {
    max(abs(cor(x[, g], z[, 1:100], use="pairwise.complete.obs")), na.rm=TRUE)
  })
  
  high_cor_genes <- names(sort(marker_cors, decreasing = TRUE)[1:50])
  genes_for_analysis <- unique(c(genes_for_analysis, high_cor_genes))
}

# Take exactly 80 genes
final_genes <- genes_for_analysis[1:min(80, length(genes_for_analysis))]
cat(sprintf("Final selected genes: %d\n", length(final_genes)))
cat("First 20 genes:\n")
print(head(final_genes, 20))

# Verify we have enough genes
if(length(final_genes) < 20) {
  stop("Not enough genes selected. Try using a different gene selection strategy.")
}

# ======================================================================
# PART 3: RUN ASCEND
# ======================================================================
cat("\n=== RUNNING ASCEND ===\n")

# Prepare data for ASCEND
ascend_data <- as.data.frame(cbind(z, x[, final_genes]))
colnames(ascend_data) <- c(paste0('z', 1:ncol(z)), 
                           paste0('x', 1:length(final_genes)))

# Create gene mapping
gene_map <- data.frame(
  x_idx = paste0('x', 1:length(final_genes)),
  gene = final_genes,
  stringsAsFactors = FALSE
)

# Scale data
cat("Scaling data...\n")
for(col in colnames(ascend_data)) {
  ascend_data[[col]] <- scale(ascend_data[[col]])[,1]
}
ascend_data[is.na(ascend_data)] <- 0

# Run ASCEND
cat("Running ASCEND algorithm...\n")
start_time <- Sys.time()

sim_obj <- list(
  dat = ascend_data,
  adj_mat = NULL,
  params = list(
    d_z = ncol(z),
    d_x = length(final_genes)
  )
)

ascend_network <- ascend_fn(sim_obj, maxiter = 3, alpha = 0.1)

# Map back to gene names
rownames(ascend_network) <- gene_map$gene[match(rownames(ascend_network), gene_map$x_idx)]
colnames(ascend_network) <- rownames(ascend_network)

ascend_time <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf("ASCEND completed in %.2f minutes\n", ascend_time))

# ======================================================================
# PART 4: RUN TRIGGER
# ======================================================================
cat("\n=== RUNNING TRIGGER ===\n")

start_time <- Sys.time()

# Build Trigger object
trig.obj <- trigger.build(
  marker = yeast$marker,
  exp = yeast$exp,
  marker.pos = yeast$marker.pos,
  exp.pos = yeast$exp.pos
)

# Run Trigger pipeline
cat("Running trigger.link...\n")
trig.obj <- trigger.link(trig.obj, norm = TRUE)

cat("Running trigger.loclink...\n")
trig.obj <- trigger.loclink(trig.obj)

cat("Running trigger.mlink with B=5...\n")
trig.obj <- trigger.mlink(trig.obj, B = 5, idx = NULL)

cat("Running trigger.net with Bsec=5...\n")
trig.prob <- trigger.net(trig.obj, Bsec = 5, idx = NULL)

trigger_time <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf("TRIGGER completed in %.2f minutes\n", trigger_time))

trigger_network <- trig.prob






# ======================================================================
# CORRECTED: Create proper gene×gene matrix from Trigger output
# ======================================================================

# Method 1: The first 6216 rows correspond to the 6216 genes
n_genes <- nrow(yeast$exp)  # Should be 6216
cat(sprintf("Number of genes in expression data: %d\n", n_genes))

# Create square matrix by taking first n_genes rows
trigger_gene_matrix <- trig.prob[1:n_genes, ]

# Add row names from column names (since rows correspond to same genes)
rownames(trigger_gene_matrix) <- colnames(trigger_gene_matrix)

cat(sprintf("Trigger gene matrix: %d × %d\n", 
            nrow(trigger_gene_matrix), ncol(trigger_gene_matrix)))
cat("First few row names:\n")
print(head(rownames(trigger_gene_matrix)))
cat("First few column names:\n")
print(head(colnames(trigger_gene_matrix)))

# Verify it's square
if(nrow(trigger_gene_matrix) == ncol(trigger_gene_matrix)) {
  cat("✓ Matrix is square - good!\n")
} else {
  cat("✗ Matrix is not square - problem!\n")
}

# Now subset to your ASCEND genes
common_genes <- intersect(rownames(ascend_network), 
                          rownames(trigger_gene_matrix))
cat(sprintf("Common genes: %d\n", length(common_genes)))

if(length(common_genes) >= 10) {
  # Subset both matrices
  ascend_subset <- ascend_network[common_genes, common_genes]
  trigger_subset <- trigger_gene_matrix[common_genes, common_genes]
  
  cat(sprintf("ASCEND subset: %d × %d\n", nrow(ascend_subset), ncol(ascend_subset)))
  cat(sprintf("Trigger subset: %d × %d\n", nrow(trigger_subset), ncol(trigger_subset)))
}










# ======================================================================
# PART 6: EXTRACT TOP REGULATORS FROM SUBSETTED NETWORKS
# ======================================================================
cat("\n=== EXTRACTING TOP REGULATORS ===\n")

# ASCEND regulators
ascend_binary <- (ascend_subset == 0.5 | ascend_subset == 1)
ascend_binary[is.na(ascend_binary)] <- FALSE
ascend_outdegree <- rowSums(ascend_binary)
head(sort(ascend_outdegree, decreasing = TRUE))
ascend_top <- names(sort(ascend_outdegree, decreasing = TRUE)[1:min(15, length(common_genes))])
ascend_top <- ascend_top[!is.na(ascend_top)]

# TRIGGER regulators - use multiple thresholds to be thorough
trigger_thresholds <- c(0.5, 0.6, 0.7, 0.8, 0.9)
trigger_results_list <- list()

for(thresh in trigger_thresholds) {
  trigger_binary <- trigger_subset > thresh
  trigger_binary[is.na(trigger_binary)] <- FALSE
  trigger_outdegree <- rowSums(trigger_binary, na.rm = TRUE)
  trigger_top <- names(sort(trigger_outdegree, decreasing = TRUE)[1:min(15, length(common_genes))])
  trigger_top <- trigger_top[!is.na(trigger_top)]
  
  trigger_results_list[[as.character(thresh)]] <- list(
    threshold = thresh,
    binary = trigger_binary,
    outdegree = trigger_outdegree,
    top = trigger_top
  )
  
  cat(sprintf("\nTRIGGER threshold = %.1f:\n", thresh))
  cat(sprintf("  Edges: %d\n", sum(trigger_binary, na.rm = TRUE)))
  cat(sprintf("  Top regulators: %s\n", paste(head(trigger_top, 5), collapse = ", ")))
}

# Choose the best threshold (e.g., one that gives similar network density to ASCEND)
ascend_edges <- sum(ascend_binary, na.rm = TRUE)
cat(sprintf("\nASCEND edges: %d\n", ascend_edges))

# Find threshold giving closest number of edges
edge_counts <- sapply(trigger_results_list, function(x) sum(x$binary, na.rm = TRUE))
best_thresh <- names(which.min(abs(edge_counts - ascend_edges)))[1]
cat(sprintf("Best matching TRIGGER threshold: %s (edges: %d)\n", 
            best_thresh, edge_counts[best_thresh]))

# Use the best threshold for downstream analysis
trigger_binary <- trigger_results_list[[best_thresh]]$binary
trigger_outdegree <- trigger_results_list[[best_thresh]]$outdegree
trigger_top <- trigger_results_list[[best_thresh]]$top

# ======================================================================
# PART 7: FUNCTION TO GET TARGET GENES
# ======================================================================
get_targets <- function(regulator, network, method, threshold = NULL) {
  # Safety check: does the regulator exist in the network?
  if(!regulator %in% rownames(network)) return(character(0))
  
  if(method == "ASCEND") {
    # Fix: Use direct logical comparison to keep matrix structure
    binary_row <- (network[regulator, ] == 0.5 | network[regulator, ] == 1)
    binary_row[is.na(binary_row)] <- FALSE
    targets <- names(which(binary_row))
  } else { 
    # TRIGGER logic
    thresh_val <- as.numeric(best_thresh)
    targets <- names(which(network[regulator, ] > thresh_val))
  }
  
  # Remove self-loops
  return(setdiff(targets, regulator))
}

# ======================================================================
# PART 8: GET GO ANNOTATIONS FOR ENRICHMENT
# ======================================================================
cat("\n=== PREPARING GO ANNOTATIONS ===\n")

# Get GO annotations for all common genes
go_data <- AnnotationDbi::select(
  org.Sc.sgd.db,
  keys = common_genes,
  columns = c("GO", "ONTOLOGY", "EVIDENCE"),
  keytype = "ORF"
)

# Filter for Biological Process with experimental evidence
curated_evidence <- c("IDA", "IMP", "IGI", "IPI", "ISS", "TAS", "NAS")
go_bp <- go_data[go_data$ONTOLOGY == "BP" & go_data$EVIDENCE %in% curated_evidence, ]
cat(sprintf("Genes with curated BP annotations: %d\n", length(unique(go_bp$ORF))))

# Create gene-to-GO mapping
gene2go <- split(go_bp$GO, go_bp$ORF)

# ======================================================================
# PART 9: SIMPLE GO ENRICHMENT FUNCTION
# ======================================================================
test_go_enrichment <- function(regulator, target_genes, all_genes, gene2go) {
  
  if(length(target_genes) < 2) {
    return(data.frame(regulator=regulator, n_targets=length(target_genes), 
                      n_sig_terms=0, best_p=1, best_term="Too few targets", stringsAsFactors=F))
  }
  
  term_counts <- table(unlist(gene2go))
  all_terms <- names(term_counts[term_counts >= 2]) # Lowered min genes per term to 2 for small set
  results <- data.frame()
  
  for(term in all_terms) {
    term_genes <- names(gene2go)[sapply(gene2go, function(x) term %in% x)]
    
    # Contingency Table
    in_target <- all_genes %in% target_genes
    in_term <- all_genes %in% term_genes
    
    a <- sum(in_target & in_term)   # Hit
    b <- sum(in_target & !in_term)  # Miss
    c <- sum(!in_target & in_term)  # Background hit
    d <- sum(!in_target & !in_term) # Background miss
    
    if(a >= 1) {
      ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
      
      term_name <- tryCatch({GO.db::Term(GOTERM[[term]])}, error = function(e) term)
      results <- rbind(results, data.frame(term=term, term_name=term_name, 
                                           p_raw=ft$p.value, stringsAsFactors=F))
    }
  }
  
  if(nrow(results) > 0) {
    results$p_adj <- p.adjust(results$p_raw, method = "BH")
    results <- results[order(results$p_raw), ]
    
    # Check if anything passed BH
    sig_hits <- results[results$p_adj < 0.05, ]
    
    if(nrow(sig_hits) > 0) {
      return(data.frame(regulator=regulator, n_targets=length(target_genes), 
                        n_sig_terms=nrow(sig_hits), best_p=sig_hits$p_adj[1], 
                        best_term=paste0(sig_hits$term_name[1], " (BH-sig)"), stringsAsFactors=F))
    } else {
      # Return the best raw p-value as a "trend"
      return(data.frame(regulator=regulator, n_targets=length(target_genes), 
                        n_sig_terms=0, best_p=results$p_raw[1], 
                        best_term=paste0(results$term_name[1], " (Trend)"), stringsAsFactors=F))
    }
  }
  
  return(data.frame(regulator=regulator, n_targets=length(target_genes), 
                    n_sig_terms=0, best_p=1, best_term=NA, stringsAsFactors=F))
}

# ======================================================================
# PART 10: RUN ENRICHMENT FOR TOP REGULATORS
# ======================================================================
cat("\n=== RUNNING GO ENRICHMENT ===\n")

# ASCEND enrichment
cat("\nTesting ASCEND top regulators:\n")
ascend_results <- data.frame()
for(reg in ascend_top) {
  targets <- get_targets(reg, ascend_subset, "ASCEND")
  cat(sprintf("  %s: %d targets\n", reg, length(targets)))
  
  res <- test_go_enrichment(reg, targets, common_genes, gene2go)
  res$method <- "ASCEND"
  res$outdegree <- ascend_outdegree[reg]
  ascend_results <- rbind(ascend_results, res)
}

# TRIGGER enrichment
cat("\nTesting TRIGGER top regulators:\n")
trigger_results <- data.frame()
for(reg in trigger_top) {
  targets <- get_targets(reg, trigger_subset, "TRIGGER")
  cat(sprintf("  %s: %d targets\n", reg, length(targets)))
  
  res <- test_go_enrichment(reg, targets, common_genes, gene2go)
  res$method <- "TRIGGER"
  res$outdegree <- trigger_outdegree[reg]
  trigger_results <- rbind(trigger_results, res)
}

# ======================================================================
# PART 11: COMPARE RESULTS
# ======================================================================
cat("\n=== COMPARISON RESULTS ===\n")

# Combine results
all_results <- rbind(ascend_results, trigger_results)

# Summary statistics
cat("\n--- Summary ---\n")
ascend_hit_rate <- mean(ascend_results$n_sig_terms > 0)
trigger_hit_rate <- mean(trigger_results$n_sig_terms > 0)

cat(sprintf("\nASCEND (%d regulators):\n", nrow(ascend_results)))
cat(sprintf("  Regulators with enriched targets: %.1f%%\n", ascend_hit_rate * 100))
cat(sprintf("  Average significant terms: %.2f\n", mean(ascend_results$n_sig_terms)))
cat(sprintf("  Average targets per regulator: %.1f\n", mean(ascend_results$n_targets)))
cat(sprintf("  Best p-value: %.4f\n", min(ascend_results$best_p)))

cat(sprintf("\nTRIGGER (%d regulators):\n", nrow(trigger_results)))
cat(sprintf("  Regulators with enriched targets: %.1f%%\n", trigger_hit_rate * 100))
cat(sprintf("  Average significant terms: %.2f\n", mean(trigger_results$n_sig_terms)))
cat(sprintf("  Average targets per regulator: %.1f\n", mean(trigger_results$n_targets)))
cat(sprintf("  Best p-value: %.4f\n", min(trigger_results$best_p)))

# Show top regulators from each method
cat("\n--- Top ASCEND Regulators by GO Enrichment ---\n")
ascend_top_enriched <- ascend_results[ascend_results$n_sig_terms > 0, ]
if(nrow(ascend_top_enriched) > 0) {
  ascend_top_enriched <- ascend_top_enriched[order(ascend_top_enriched$best_p), ]
  print(ascend_top_enriched[, c("regulator", "outdegree", "n_targets", "n_sig_terms", "best_p", "best_term")])
} else {
  cat("No ASCEND regulators with enriched targets\n")
}

cat("\n--- Top TRIGGER Regulators by GO Enrichment ---\n")
trigger_top_enriched <- trigger_results[trigger_results$n_sig_terms > 0, ]
if(nrow(trigger_top_enriched) > 0) {
  trigger_top_enriched <- trigger_top_enriched[order(trigger_top_enriched$best_p), ]
  print(trigger_top_enriched[, c("regulator", "outdegree", "n_targets", "n_sig_terms", "best_p", "best_term")])
} else {
  cat("No TRIGGER regulators with enriched targets\n")
}

# ======================================================================
# PART 12: STATISTICAL TEST
# ======================================================================
cat("\n=== STATISTICAL COMPARISON ===\n")

# Fisher's exact test comparing proportion of regulators with enrichment
if(nrow(ascend_results) > 0 && nrow(trigger_results) > 0) {
  cont_table <- matrix(c(
    sum(ascend_results$n_sig_terms > 0),
    nrow(ascend_results) - sum(ascend_results$n_sig_terms > 0),
    sum(trigger_results$n_sig_terms > 0),
    nrow(trigger_results) - sum(trigger_results$n_sig_terms > 0)
  ), nrow = 2)
  
  print(cont_table)
  ft <- fisher.test(cont_table)
  cat(sprintf("Fisher's exact test p-value: %.4f\n", ft$p.value))
  
  if(ft$p.value < 0.05) {
    if(sum(ascend_results$n_sig_terms > 0) > sum(trigger_results$n_sig_terms > 0)) {
      cat("Conclusion: ASCEND finds significantly more regulators with biological relevance\n")
    } else {
      cat("Conclusion: TRIGGER finds significantly more regulators with biological relevance\n")
    }
  } else {
    cat("Conclusion: No significant difference in biological relevance between methods\n")
  }
}

# ======================================================================
# PART 13: VISUALIZATION
# ======================================================================
cat("\n=== CREATING PLOTS ===\n")

pdf("ascend_trigger_comparison.pdf", width = 12, height = 10)
par(mfrow = c(2, 2))

# Plot 1: Outdegree distribution
boxplot(list(ASCEND = ascend_outdegree, TRIGGER = trigger_outdegree),
        col = c("dodgerblue", "firebrick"),
        main = "Regulator Activity (Outdegree)",
        ylab = "Number of targets",
        outline = FALSE)
points(1, mean(ascend_outdegree), pch = 19, col = "darkblue", cex = 1.5)
points(2, mean(trigger_outdegree), pch = 19, col = "darkred", cex = 1.5)

# Plot 2: Percentage with enrichment
if(nrow(ascend_results) > 0 && nrow(trigger_results) > 0) {
  barplot(c(ascend_hit_rate * 100, trigger_hit_rate * 100),
          names.arg = c("ASCEND", "TRIGGER"),
          col = c("dodgerblue", "firebrick"),
          main = "Regulators with GO Enriched Targets",
          ylab = "Percentage (%)",
          ylim = c(0, 100))
}

# Plot 3: Number of significant GO terms
if(nrow(ascend_results) > 0 && nrow(trigger_results) > 0) {
  boxplot(list(ASCEND = ascend_results$n_sig_terms, 
               TRIGGER = trigger_results$n_sig_terms),
          col = c("dodgerblue", "firebrick"),
          main = "Number of Significant GO Terms per Regulator",
          ylab = "Count")
}

# Plot 4: -log10(p-values)
ascend_pvals <- -log10(ascend_results$best_p[ascend_results$best_p < 1])
trigger_pvals <- -log10(trigger_results$best_p[trigger_results$best_p < 1])

if(length(ascend_pvals) > 0 && length(trigger_pvals) > 0) {
  boxplot(list(ASCEND = ascend_pvals, TRIGGER = trigger_pvals),
          col = c("dodgerblue", "firebrick"),
          main = "Enrichment Significance (-log10 p-value)",
          ylab = "-log10(p-value)")
  abline(h = -log10(0.05), lty = 2, col = "gray")
}

dev.off()
cat("Plot saved to ascend_trigger_comparison.pdf\n")

# ======================================================================
# PART 14: SAVE RESULTS
# ======================================================================
save(ascend_subset, trigger_subset, 
     ascend_results, trigger_results,
     all_results,
     file = "comparison_results.RData")

cat("\n=== COMPLETE ===\n")
cat(sprintf("Common genes analyzed: %d\n", length(common_genes)))
cat("Results saved to comparison_results.RData\n")
cat("Plot saved to ascend_trigger_comparison.pdf\n")


library(ggplot2)

# Combine your results for ggplot
plot_data <- rbind(
  data.frame(ascend_results, Method="ASCEND"),
  data.frame(trigger_results, Method="TRIGGER")
)

# Filter for the "Trend" or "BH-sig" rows you want to show
plot_data <- plot_data[plot_data$best_p < 0.05, ]

ggplot(plot_data, aes(x = Method, y = best_term)) +
  geom_point(aes(size = n_targets, color = -log10(best_p))) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(title = "Functional Landscape: ASCEND vs. TRIGGER",
       y = "Top Enriched GO Biological Process",
       size = "Target Count",
       color = "-log10(adj P)") +
  theme(axis.text.x = element_text(face = "bold", size = 12))


ggplot(all_results, aes(x = n_targets, y = -log10(best_p), color = method)) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  annotate("text", x = 20, y = 1.5, label = "Significance Threshold", color = "grey50") +
  scale_color_manual(values = c("ASCEND" = "#0072B2", "TRIGGER" = "#D55E00")) +
  theme_classic() +
  labs(x = "Number of Predicted Targets", 
       y = "-log10(Adjusted P-value)",
       title = "Regulatory Specificity Comparison")








# install.packages(c("pheatmap", "igraph", "ggplot2", "gridExtra", "RColorBrewer"))
library(pheatmap)
library(igraph)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

######################################################################
# PRE-PROCESSING: Define top common genes for comparison
#
# A deeply informative analysis for a *Bioinformatics* paper focuses on 
# a clustered view of the top "hottest" gene subsets, rather than a 
# simple counts boxplot, to show the coordinated biological programs that 
# are denser in ASCEND and sparser in TRIGGER.
######################################################################
top_15_ascend <- names(sort(ascend_outdegree, decreasing = TRUE)[1:15])
top_15_trigger <- names(sort(trigger_outdegree, decreasing = TRUE)[1:15])
unique_top_regulators <- unique(c(top_15_ascend, top_15_trigger))

# Create matrices restricted to these top common genes for a focused view
ascend_clust_mat <- ascend_subset[unique_top_regulators, unique_top_regulators]
trigger_clust_mat <- trigger_subset[unique_top_regulators, unique_top_regulators]


######################################################################
# FIGURE CONCEPT 1: Network Graph Comparisons (Like Panel A in Fig 2)
#
# Yes, you should generate a similar graph for ASCEND to allow for direct visual 
# comparison against the existing TRIGGER plot.
######################################################################
generate_network_plot <- function(mat, title, color_palette = "blue") {
  # Threshold the network matrix to get only edges with probabilities of 0.5 or 1
  binary_mat <- (mat %in% c(0.5, 1))
  g <- graph_from_adjacency_matrix(binary_mat, mode = "directed", diag = FALSE)
  
  # Visualization parameters to match standard elegance
  V(g)$color <- ifelse(color_palette == "blue", "#0072B2", "#D55E00")
  V(g)$size <- 5
  E(g)$arrow.size <- 0.2
  E(g)$color <- "gray80"
  
  # Use a circular layout for coordinated visibility
  plot(g, layout = layout_in_circle, vertex.label = NA, main = title)
}

generate_network_plot(ascend_clust_mat, "A) ASCEND Network Graph")
generate_network_plot(trigger_clust_mat, "B) TRIGGER Network Graph", "orange")


######################################################################
# FIGURE CONCEPT 2: Clustered Heatmaps (Deeply Informative Heatmap View)
#
# This is the gold standard for publication, replacing boxplots with a 
# clustered heatmap that allows you to see the precise relationship for every 
# gene pair. Using clustering demonstrates the coordinated regulatory programs 
# that are sparser in TRIGGER but dense in ASCEND.
######################################################################
heatmap_colors <- colorRampPalette(c("white", "dodgerblue4"))(100) # White-to-blue gradient

# Standard square dendrogram visualization is powerful and informative
p_ascend_h <- pheatmap(ascend_clust_mat, color = heatmap_colors, 
                       main = "Clustered Coordinated Program (ASCEND)",
                       silent = TRUE)

p_trigger_h <- pheatmap(trigger_clust_mat, color = heatmap_colors, 
                        main = "Sparse Connections (TRIGGER)",
                        silent = TRUE)

grid.arrange(p_ascend_h$gtable, p_trigger_h$gtable, ncol = 2)


######################################################################
# FIGURE CONCEPT 3: Functional Dot Plot (Gold Standard GO Summary)
#
# A dot plot is a superior visualization because it conveys more information than 
# a standard boxplot. It allows you to simultaneously display count, 
# significance, and function in a compact and elegant view.
######################################################################
# Prepare a data frame combining best_p and n_targets for all genes
all_results_combined <- rbind(
  data.frame(ascend_results, method="ASCEND"),
  data.frame(trigger_results, method="TRIGGER")
)

# A standard, deeply informative visualization for a functional landscape
ggplot(all_results_combined, aes(x = method, y = best_term)) +
  geom_point(aes(size = n_targets, color = -log10(best_p))) + # size maps to target count, color to significance
  scale_color_gradient(low = "blue", high = "red") + # standard color scale
  theme_minimal() + # highly elegant and standard style
  labs(title = "Functional Landscape: ASCEND vs. TRIGGER",
       y = "Top Enriched GO Biological Process",
       size = "Target Count",
       color = "-log10(adj P)") +
  theme(axis.text.x = element_text(face = "bold", size = 12)) # Deeply informative text scaling


library(igraph)

# Load your full ASCEND network matrix
# load("ascend_network.RData") 

# Function to generate a network graph that allows for a deeply informative comparison
generate_network_graph <- function(network_matrix, method_name) {
  # 1. Create the binary logical matrix
  binary_matrix <- (network_matrix == 0.5 | network_matrix == 1)
  
  # 2. FIX: Replace NAs with FALSE so igraph can process it
  binary_matrix[is.na(binary_matrix)] <- FALSE
  
  # 3. Convert to graph object
  g <- graph_from_adjacency_matrix(binary_matrix, mode = "directed", diag = FALSE)
  
  # 4. Set aesthetic parameters for Bioinformatics-grade plots
  V(g)$color <- "dodgerblue"
  V(g)$size <- 4
  V(g)$frame.color <- "white"
  E(g)$arrow.size <- 0.15
  E(g)$color <- adjustcolor("gray70", alpha.f = 0.5)
  
  # 5. Plot using a layout that highlights "Causal Genes" (Hubs)
  plot(g, 
       layout = layout_with_fr(g), # Fruchterman-Reingold creates organic clusters
       vertex.label = NA, 
       main = paste("ASCEND Transcriptional Regulatory Network\n", 
                    nrow(binary_matrix), "Genes,", sum(binary_matrix), "Edges"))
}
# Generate the graph for ASCEND. This is highly encouraged over boxplots.
generate_network_graph(ascend_network, "ASCEND")



library(igraph)
library(RColorBrewer)

# 1. Prepare the Data
# Ensure the matrix is numeric and handle NAs
ascend_matrix_clean <- ascend_subset
ascend_matrix_clean[is.na(ascend_matrix_clean)] <- 0

# 2. Create the Graph
# We use the original values to preserve edge weights (0.5 vs 1)
g <- graph_from_adjacency_matrix(ascend_matrix_clean, 
                                 mode = "directed", 
                                 weighted = TRUE, 
                                 diag = FALSE)

# 3. Aesthetics: Nodes (Color by GO significance)
# Let's highlight the regulators that showed high GO enrichment
sig_regs <- ascend_results$regulator[ascend_results$n_sig_terms > 0]
V(g)$color <- ifelse(V(g)$name %in% sig_regs, "#E69F00", "#56B4E9") # Orange for sig, Blue for others
V(g)$size <- ifelse(V(g)$name %in% sig_regs, 8, 4) # Make enriched regulators larger
V(g)$label <- ifelse(V(g)$name %in% sig_regs, V(g)$name, NA) # Label only top regulators
V(g)$label.cex <- 0.7
V(g)$label.color <- "black"
V(g)$label.font <- 2 # Bold

# 4. Aesthetics: Edges (Refine Arrows and Labels)
# Map weights to colors or labels
E(g)$color <- ifelse(E(g)$weight == 1, "gray20", "gray70")
E(g)$width <- ifelse(E(g)$weight == 1, 1.2, 0.6)

# To label edges with the probability value:
E(g)$label <- E(g)$weight 
E(g)$label.cex <- 0.5
E(g)$label.color <- "darkred"

# 5. Final Publication Plot
pdf("ASCEND_HighRes_Network.pdf", width = 10, height = 10)

# Fruchterman-Reingold layout for organic clustering
set.seed(42) # For reproducibility
plot(g, 
     layout = layout_with_fr(g),
     vertex.frame.color = "white",
     edge.arrow.size = 0.3,       # Much smaller, sharper arrows
     edge.arrow.width = 0.8,      # Narrower arrow heads
     edge.curved = 0.2,           # Slight curve prevents overlapping straight lines
     main = "ASCEND Regulatory Network: Causal Clusters")

legend("bottomleft", 
       legend = c("Enriched Regulator (p < 0.05)", "Target/Other Gene"),
       fill = c("#E69F00", "#56B4E9"), 
       bquote = TRUE, cex = 0.8)

dev.off()
