library(readr)
rna_counts <- read_csv("GSE173380_RNAseq_counts.csv")
head(rna_counts)

colnames(rna_counts)


library(dplyr)
library(readr)
library(stringr)
library(purrr)

meth_path <- "gse173380_meth/"
meth_files <- list.files(meth_path, pattern = "*.txt", full.names = TRUE)

# Read all methylation files and tag each with its sample name

# selecting the right header patterns ie BN_NP1, CD_P2, etc.
pattern <- "(?:BN|CD)_(?:NP|P)\\d+"

# quick test of regex on basenames
basename(meth_files)
str_extract(basename(meth_files), pattern)

# Read all methylation files and tag each with its sample name (robust)
meth_all <- map_dfr(meth_files, function(f) {
  b <- basename(f)
  sample_id <- str_extract(b, pattern)
  if (is.na(sample_id)) {
    warning("Could not extract sample id from filename: ", b)
  }
  df <- read_tsv(f, show_col_types = FALSE)
  df$sample <- sample_id
  df$source_file <- b
  return(df)
})















########################################################################
## MAIN CODE FOR PROCESSING METHYLATION DATA FILES ########################
########################################################################

#### START: Full ASCEND-ready preprocessing pipeline (R) ####

## ---------- PACKAGES (install if needed) ----------
# Uncomment installs if you need them:
# if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
# BiocManager::install(c("biomaRt","GenomicRanges","IRanges","S4Vectors","GenomicFeatures"))
# install.packages(c("readr","dplyr","tidyr","stringr","purrr"))

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(biomaRt)
library(GenomicRanges)
library(IRanges)

## ---------- USER PARAMETERS (adjust as needed) ----------
meth_dir <- "gse173380_meth"                # directory with GSM*.txt methylation files
rna_counts_file <- "GSE173380_RNAseq_counts.csv"  # RNA counts CSV (first col = Ensembl gene id)
prom_upstream <- 2000       # promoter upstream (bp)
prom_downstream <- 500      # promoter downstream (bp)
cov_thresh <- 5             # minimum coverage per CpG to keep
min_samples_with_meth <- 3  # minimum number of samples with methylation required to keep gene
var_threshold_expr <- 0.2   # min SD on expression (z-scored thresholding done later)
var_threshold_meth <- 0.02  # min SD on methylation
pseudo <- 1e-6              # small pseudo-count to avoid logit infinities

## ---------- 1) Read RNA counts ----------
rna_counts <- read_csv(rna_counts_file, show_col_types = FALSE)
# Ensure first column is gene id
names(rna_counts)[1] <- "ensembl_gene_id"
# Strip version numbers if present
rna_counts <- rna_counts %>%
  mutate(ensembl_gene_id = sub("\\.\\d+$", "", ensembl_gene_id))

message("RNA counts: ", nrow(rna_counts), " genes, ", ncol(rna_counts)-1, " samples.")

## ---------- 2) Read methylation (RRBS) files and combine ----------
meth_files <- list.files(meth_dir, pattern = "\\.txt$", full.names = TRUE)
if(length(meth_files) == 0) stop("No methylation .txt files found in ", meth_dir)

# helper to extract sample id from filename (expects BN_/CD_ pattern)
extract_sample <- function(fname) {
  b <- basename(fname)
  s <- str_extract(b, "(?:BN|CD)_(?:NP|P)\\d+")
  if(is.na(s)) {
    # fallback: split by '_' and take last two parts
    parts <- unlist(str_split(b, "_"))
    if(length(parts) >= 3) {
      s <- paste0(parts[length(parts)-1], "_", str_remove(parts[length(parts)], "\\.txt$"))
    }
  }
  return(s)
}

# read and tag
read_one_meth <- function(f) {
  df <- read_tsv(f, show_col_types = FALSE, col_types = cols(.default = "c"))
  # expected columns: chrBase chr base strand coverage freqC freqT
  # coerce typed columns
  if(!("base" %in% names(df)) || !("chr" %in% names(df))) {
    stop("Unexpected columns in methylation file ", f, ". Inspect file headers.")
  }
  df <- df %>%
    mutate(
      chr = as.character(chr),
      base = as.integer(base),
      coverage = as.numeric(coverage),
      freqC = as.numeric(freqC),
      sample = extract_sample(f),
      source_file = basename(f)
    ) %>%
    # --- ADD THIS FIX HERE ---
    mutate(
      strand = case_match(
        strand,
        "1" ~ "+",
        "+" ~ "+",
        "-1" ~ "-",
        "-" ~ "-",
        .default = "*" # Set anything unexpected to * (unknown/ambiguous)
      )
    ) %>%
    dplyr::select(chr, base, strand, coverage, freqC, sample, source_file)
  return(df)
}

message("Reading methylation files (this may take a few seconds)...")
meth_list <- map(meth_files, read_one_meth)
meth_all <- bind_rows(meth_list)

# quick check
table(meth_all$sample, useNA = "ifany")

## ---------- 3) Filter low-coverage CpGs and compute methylation fraction ----------
meth_all <- meth_all %>%
  filter(!is.na(coverage) & coverage >= cov_thresh) %>%
  mutate(methylation = pmin(pmax(freqC / 100, 0), 1)) # fraction 0..1

message("After coverage filter: CpG rows = ", nrow(meth_all),
        " ; distinct samples = ", length(unique(meth_all$sample)))

## ---------- 4) Build CpG GRanges (note chr prefix) ----------
# Add 'chr' prefix if not present (common in these files: '10' etc)
if(!grepl("^chr", meth_all$chr[1], ignore.case = TRUE)) {
  meth_all <- meth_all %>% mutate(chr = paste0("chr", chr))
}
meth_gr <- GRanges(seqnames = meth_all$chr,
                   ranges = IRanges(start = meth_all$base, end = meth_all$base),
                   strand = meth_all$strand,
                   methylation = meth_all$methylation,
                   coverage = meth_all$coverage,
                   sample = meth_all$sample)

## ---------- 5) Obtain gene coordinates from Ensembl via biomaRt ----------
message("Querying Ensembl (biomaRt) for Rattus norvegicus gene coordinates...")
ensembl <- useEnsembl(biomart = "genes", dataset = "rnorvegicus_gene_ensembl", mirror = "useast")

gene_ids <- unique(rna_counts$ensembl_gene_id)
# Chunking if too large might be needed; we'll query them directly
bm <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position",
                           "end_position", "strand"),
            filters = "ensembl_gene_id",
            values = gene_ids,
            mart = ensembl)

if(nrow(bm) == 0) stop("biomaRt returned no gene coordinates - check network or gene ID types.")

# standardize chromosome names to 'chrN'
bm <- bm %>%
  mutate(chromosome_name = as.character(chromosome_name)) %>%
  mutate(seqname = ifelse(grepl("^chr", chromosome_name, ignore.case = TRUE),
                          chromosome_name, paste0("chr", chromosome_name)))

# compute promoter coordinates by strand (TSS +/-)
bm <- bm %>%
  rowwise() %>%
  mutate(prom_start = ifelse(strand == 1,
                             pmax(1, start_position - prom_upstream),
                             pmax(1, end_position - prom_downstream)),
         prom_end = ifelse(strand == 1,
                           start_position + prom_downstream,
                           end_position + prom_upstream)) %>%
  ungroup()

prom_gr <- GRanges(seqnames = bm$seqname,
                   ranges = IRanges(start = bm$prom_start, end = bm$prom_end),
                   strand = ifelse(bm$strand == 1, "+", "-"),
                   ensembl_gene_id = bm$ensembl_gene_id)

message("Promoters prepared: ", length(prom_gr), " entries.")

## ---------- 6) Overlap CpGs -> promoters and aggregate (coverage-weighted) ----------
hits <- findOverlaps(prom_gr, meth_gr, ignore.strand = FALSE)
if(length(hits) == 0) stop("No overlaps found between CpGs and promoters. Check genome build / chr names.")

ov_df <- data.frame(prom_idx = queryHits(hits),
                    cpg_idx = subjectHits(hits),
                    gene = prom_gr$ensembl_gene_id[queryHits(hits)],
                    methylation = mcols(meth_gr)$methylation[subjectHits(hits)],
                    coverage = mcols(meth_gr)$coverage[subjectHits(hits)],
                    sample = mcols(meth_gr)$sample[subjectHits(hits)],
                    stringsAsFactors = FALSE)

# compute weighted methylation per gene per sample
agg <- ov_df %>%
  group_by(gene, sample) %>%
  summarize(n_cpg = n(),
            total_cov = sum(coverage, na.rm = TRUE),
            meth_weighted = ifelse(total_cov > 0,
                                   sum(methylation * coverage, na.rm = TRUE) / total_cov,
                                   mean(methylation, na.rm = TRUE)),
            meth_mean = mean(methylation, na.rm = TRUE),
            .groups = "drop")

# wide matrix: gene x sample
meth_wide <- agg %>%
  dplyr::select(gene, sample, meth_weighted) %>%
  pivot_wider(names_from = sample, values_from = meth_weighted)

# ensure row order matches RNA gene list; strip gene version if needed
rownames_meth <- as.character(meth_wide$gene)
meth_wide_df <- as.data.frame(meth_wide)
rownames(meth_wide_df) <- meth_wide_df$gene
meth_wide_df$gene <- NULL

# reorder rows to RNA gene order (if gene missing, row will be NA)
rna_gene_order <- rna_counts$ensembl_gene_id
rna_gene_order_clean <- sub("\\.\\d+$", "", rna_gene_order)
rownames(meth_wide_df) <- sub("\\.\\d+$", "", rownames(meth_wide_df))
meth_final <- meth_wide_df[rna_gene_order_clean, , drop = FALSE]
# add ensembl id column back (with original RNA IDs)
meth_final_df <- data.frame(ensembl_gene_id = rna_gene_order, meth_final, check.names = FALSE)

## ---------- 7) Logit transform methylation (handle 0/1 with pseudo) ----------
logit_fun <- function(p) log((p + pseudo) / (1 - p + pseudo))
meth_numeric <- as.matrix(meth_final_df[,-1])
# If matrix is entirely NA for some columns, leave as NA
meth_logit <- apply(meth_numeric, 2, function(col) {
  if(all(is.na(col))) return(rep(NA_real_, length(col)))
  return(logit_fun(col))
})
# convert back to data.frame with gene ids
meth_logit_df <- data.frame(ensembl_gene_id = meth_final_df$ensembl_gene_id,
                            as.data.frame(meth_logit, check.names = FALSE),
                            check.names = FALSE, stringsAsFactors = FALSE)

## ---------- 8) Normalize RNA: CPM -> log2CPM ----------
counts_mat <- as.matrix(rna_counts[ , -1])
rownames(counts_mat) <- rna_counts$ensembl_gene_id
libsize <- colSums(counts_mat, na.rm = TRUE)
cpm <- sweep(counts_mat, 2, libsize/1e6, FUN = "/")
log2cpm <- log2(cpm + 1)
log2cpm_df <- data.frame(ensembl_gene_id = rownames(log2cpm), as.data.frame(log2cpm, check.names = FALSE), check.names = FALSE)

## ---------- 9) Z-score (standardize) both matrices across samples (columns) ----------
# helper that z-scores columns (center=TRUE, scale=TRUE). Keep gene column.
zscale_df <- function(df){
  idcol <- df[[1]]
  mat <- as.matrix(df[,-1])
  # compute column means/sds or use scale
  mat_z <- apply(mat, 2, function(x) {
    # if all NA, keep NA
    if(all(is.na(x))) return(x)
    # else z-score (center & scale) but preserve NA
    mu <- mean(x, na.rm = TRUE)
    s  <- sd(x, na.rm = TRUE)
    if(is.na(s) || s == 0) return((x - mu))
    (x - mu) / s
  })
  out <- data.frame(ensembl_gene_id = idcol, as.data.frame(mat_z, check.names = FALSE), check.names = FALSE, stringsAsFactors = FALSE)
  return(out)
}

expr_z <- zscale_df(log2cpm_df)
meth_z  <- zscale_df(meth_logit_df)

## ---------- 10) Filter genes by methylation presence & variance ----------
# Keep genes that have methylation in at least min_samples_with_meth samples
meth_nonNA_counts <- rowSums(!is.na(as.matrix(meth_z[,-1])))
keep_by_presence <- meth_nonNA_counts >= min_samples_with_meth

# Compute SD across samples for each gene (on z-scored matrices)
expr_sd <- apply(as.matrix(expr_z[,-1]), 1, sd, na.rm = TRUE)
meth_sd <- apply(as.matrix(meth_z[,-1]), 1, sd, na.rm = TRUE)

keep_by_variance <- (expr_sd >= var_threshold_expr) & (meth_sd >= var_threshold_meth)

keep_genes <- which(keep_by_presence & keep_by_variance)
kept_gene_ids <- expr_z$ensembl_gene_id[keep_genes]

message(length(kept_gene_ids), " genes kept after presence & variance filtering.")

expr_final <- expr_z %>% filter(ensembl_gene_id %in% kept_gene_ids)
meth_final2 <- meth_z %>%   filter(ensembl_gene_id %in% kept_gene_ids)

## ---------- 11) Align sample order in both matrices ----------
# identify common sample columns (should be same)
expr_samples <- colnames(expr_final)[-1]
meth_samples <- colnames(meth_final2)[-1]
common_samples <- intersect(expr_samples, meth_samples)
if(length(common_samples) == 0) stop("No common sample names between RNA and methylation matrices. Check sample naming.")

# Subset & reorder columns consistently
expr_ready <- expr_final %>%
  dplyr::select(ensembl_gene_id, all_of(common_samples))
meth_ready <- meth_final2 %>%
  dplyr::select(ensembl_gene_id, all_of(common_samples))

# Final check: identical sample column order
if(!identical(colnames(expr_ready)[-1], colnames(meth_ready)[-1])) {
  stop("Sample column orders differ after alignment. Inspect column names.")
}

## ---------- 12) Save outputs for ASCEND ----------
write_csv(expr_ready, "EXPRESSION_log2CPM_zscored.csv")
write_csv(meth_ready, "promoter_METHYLATION_logit_zscored.csv")

# Also produce merged wide table (interleaved: expr cols then meth cols)
merged <- inner_join(expr_ready, meth_ready, by = "ensembl_gene_id", suffix = c("_expr","_meth"))
write_csv(merged, "MERGED_expression_and_methylation.csv")

message("Wrote 3 files:\n - ASCEND_expression_log2CPM_zscored.csv\n - ASCEND_promoter_methylation_logit_zscored.csv\n - ASCEND_merged_expression_and_methylation.csv")




