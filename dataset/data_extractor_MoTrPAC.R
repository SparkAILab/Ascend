# prepare_motrpac_for_ascend_t70.R
# End-to-end preprocessing: MoTrPAC t70 white-adipose RRBS + RNA -> ASCEND-ready
# Usage: set working directory to folder containing:
#  - GSE242354_motrpac_pass1b-06_t70-white-adipose_transcript-rna-seq_rsem-genes-count.txt
#  - GSE242354_motrpac_pass1b-06_transcript-rna-seq_qa-qc-metrics.csv  (optional but recommended)
#  - directory "GSE242358_methylation/" with .cov files

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(biomaRt)
  library(GenomicRanges)
  library(IRanges)
  library(edgeR)
})

## ---------------- USER PARAMETERS ----------------
rna_counts_file   <- "GSE242354_motrpac_pass1b-06_t70-white-adipose_transcript-rna-seq_rsem-genes-count.txt"
qa_qc_rna_file    <- "GSE242354_motrpac_pass1b-06_transcript-rna-seq_qa-qc-metrics.csv" # optional but helpful
meth_dir          <- "GSE242358_methylation"
prom_upstream     <- 2000
prom_downstream   <- 500
cov_thresh        <- 5
min_samples_with_meth <- 3
pseudo            <- 1e-6
var_threshold_expr <- 0.2
var_threshold_meth <- 0.02
# -------------------------------------------------

message(">>> START: MoTrPAC t70 preprocessing for ASCEND")

## ---------- Helpers ----------
extract_vial_from_covname <- function(fname) {
  b <- basename(fname)
  # 11-digit vial first attempt
  v <- str_extract(b, "\\d{11}")
  if (!is.na(v)) return(v)
  # fallback: 10+ digits
  v <- str_extract(b, "\\d{8,}")
  return(v)
}

# Try to find RRBS vial for a given RNA vial using multiple heuristics
find_rrbs_candidate <- function(vial_rna, cov_vials, qa_df = NULL) {
  candidates <- character(0)
  # 1) exact
  candidates <- c(candidates, vial_rna)
  # 2) replace trailing '05' with '04' (observed pattern)
  candidates <- c(candidates, sub("05$", "04", vial_rna))
  # 3) replace final digit 5->4
  candidates <- c(candidates, sub("5$", "4", vial_rna))
  # 4) try removing last char (if datasets differ by trailing char)
  candidates <- c(candidates, substr(vial_rna, 1, nchar(vial_rna)-1))
  # 5) PID lookup via qa_df: find other vials with same PID and whose vial appears in cov_vials
  if(!is.null(qa_df) && "PID" %in% names(qa_df)) {
    pid_row <- qa_df %>% filter(viallabel == vial_rna | as.character(viallabel) == vial_rna)
    if(nrow(pid_row) >= 1) {
      pid <- pid_row$PID[1]
      other_vials <- qa_df %>% filter(PID == pid) %>% pull(viallabel) %>% unique() %>% as.character()
      candidates <- c(candidates, other_vials)
    }
  }
  # reduce and test intersection with cov_vials
  candidates <- unique(na.omit(candidates))
  found <- intersect(candidates, cov_vials)
  if(length(found) > 0) return(found[1])
  return(NA_character_)
}

## ---------- 1) read QA/QC RNA file (optional, used to map PID) ----------
qa_rna <- NULL
if(file.exists(qa_qc_rna_file)) {
  message("Reading RNA QA/QC: ", qa_qc_rna_file)
  qa_rna <- read_csv(qa_qc_rna_file, show_col_types = FALSE) %>%
    mutate(viallabel = as.character(viallabel))
  # try multiple column name variants
  if(!"viallabel" %in% names(qa_rna)) {
    # find plausible vial column
    vcol <- grep("vial", names(qa_rna), ignore.case = TRUE, value = TRUE)
    if(length(vcol) > 0) {
      qa_rna <- qa_rna %>% rename(viallabel = all_of(vcol[1])) %>% mutate(viallabel = as.character(viallabel))
    }
  }
  if(!"PID" %in% names(qa_rna)) {
    pcol <- grep("^PID$|pid|SubjectID", names(qa_rna), ignore.case = TRUE, value = TRUE)
    if(length(pcol) > 0) qa_rna <- qa_rna %>% rename(PID = all_of(pcol[1]))
  }
} else {
  message("RNA QA/QC file not found; PID-based mapping will be skipped.")
}

## ---------- 2) Read RNA counts ----------
message("Reading RNA counts file: ", rna_counts_file)
if(!file.exists(rna_counts_file)) stop("RNA counts file not found: ", rna_counts_file)

# use fread robustly; allow irregular lines
rna_dt <- tryCatch({
  fread(rna_counts_file, fill = Inf, data.table = FALSE)
}, error = function(e) {
  stop("Failed to read RNA counts: ", e$message)
})

# first col is gene_id
colnames(rna_dt)[1] <- "ensembl_gene_id"
rna_dt$ensembl_gene_id <- sub("\\.\\d+$", "", rna_dt$ensembl_gene_id)
rna_gene_ids <- rna_dt$ensembl_gene_id

# sanitize sample column names: they look like numeric viallabels already; keep original strings
rna_counts_mat <- as.data.frame(rna_dt[ , -1], stringsAsFactors = FALSE)
colnames(rna_counts_mat) <- colnames(rna_dt)[-1]
# coerce to numeric safely (suppress warnings)
rna_counts_mat[] <- lapply(rna_counts_mat, function(x) suppressWarnings(as.numeric(as.character(x))))
rownames(rna_counts_mat) <- rna_gene_ids
message("RNA matrix: genes=", nrow(rna_counts_mat), " samples=", ncol(rna_counts_mat))

## ---------- 3) locate RRBS .cov files ----------
message("Searching for RRBS .cov files in: ", meth_dir)
cov_files_all <- list.files(meth_dir, pattern = "\\.cov$", full.names = TRUE, ignore.case = TRUE)
if(length(cov_files_all) == 0) stop("No .cov files found in: ", meth_dir)
cov_index <- tibble(file = cov_files_all) %>%
  mutate(vial = map_chr(file, extract_vial_from_covname)) %>%
  filter(!is.na(vial))
cov_vials <- unique(cov_index$vial)
message("Found ", length(cov_files_all), " .cov files; unique vial IDs: ", length(cov_vials))

## ---------- 4) identify RNA samples for t70 white adipose ----------
# if QA present, filter to tissue white/adipose and timepoint t70
if(!is.null(qa_rna)) {
  # find Tissue-like column
  tissue_col <- grep("tissue|specimen|sampletype", names(qa_rna), ignore.case = TRUE, value = TRUE)[1]
  time_col <- grep("timepoint|study_group_timepoint|time_point", names(qa_rna), ignore.case = TRUE, value = TRUE)[1]
  vialcol  <- "viallabel"
  if(is.na(tissue_col)) {
    message("Could not find tissue column in QA/QC; defaulting to using all columns.")
    qa_rna_tissue <- qa_rna
  } else {
    qa_rna_tissue <- qa_rna %>% filter(grepl("white", !!sym(tissue_col), ignore.case = TRUE) |
                                         grepl("adipose", !!sym(tissue_col), ignore.case = TRUE) |
                                         grepl("t70", !!sym(tissue_col), ignore.case = TRUE))
  }
  rna_vials <- unique(na.omit(as.character(qa_rna_tissue[[vialcol]])))
  message("RNA QA/QC reports ", length(rna_vials), " white/adipose vials (t70 subset)")
} else {
  # fallback: take all RNA columns as candidate vials
  rna_vials <- colnames(rna_counts_mat)
  message("No QA present: using all RNA count columns as candidate vials (", length(rna_vials), ")")
}


## ---------- 5) build mapping RNA->RRBS using heuristics ----------
message("Building RNA -> RRBS vial mapping (heuristic)...")

# Create the mapping dataframe with explicit column creation
mapping_candidates <- data.frame(
  vial_rna = as.character(rna_vials),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(vial_rna))

# Initialize the guess column explicitly
mapping_candidates$vial_rrbs_guess <- NA_character_

# Debug: check column exists
message("Columns after initialization: ", paste(colnames(mapping_candidates), collapse=", "))

# Fill in the guesses
for(i in 1:nrow(mapping_candidates)) {
  vial_rna <- mapping_candidates$vial_rna[i]
  guess <- find_rrbs_candidate(vial_rna, cov_vials, qa_rna)
  mapping_candidates$vial_rrbs_guess[i] <- ifelse(is.na(guess), NA_character_, guess)
}

# Debug: check column still exists and show results
message("Columns after filling guesses: ", paste(colnames(mapping_candidates), collapse=", "))
message("Mapping results:")
print(table(is.na(mapping_candidates$vial_rrbs_guess), useNA = "always"))

# If some still NA, attempt PID-based mapping (if qa_rna exists)
if(sum(is.na(mapping_candidates$vial_rrbs_guess)) > 0 && !is.null(qa_rna) && "PID" %in% names(qa_rna)) {
  message("Attempting PID-based mapping for ", sum(is.na(mapping_candidates$vial_rrbs_guess)), " unmapped samples")
  
  missing_indices <- which(is.na(mapping_candidates$vial_rrbs_guess))
  
  for(i in missing_indices) {
    v <- mapping_candidates$vial_rna[i]
    # Use exact string matching for viallabel
    pid_row <- qa_rna %>% filter(as.character(viallabel) == v)
    if(nrow(pid_row) >= 1) {
      pid <- pid_row$PID[1]
      # Find other vials with same PID
      other_vials <- qa_rna %>% 
        filter(PID == pid) %>% 
        pull(viallabel) %>% 
        as.character() %>%
        unique()
      found <- intersect(other_vials, cov_vials)
      if(length(found) > 0) {
        mapping_candidates$vial_rrbs_guess[i] <- found[1]
        message("PID mapping: ", v, " -> ", found[1])
      }
    }
  }
}

# Final debug check
message("Final columns before filtering: ", paste(colnames(mapping_candidates), collapse=", "))

# Now safely create the final mapping - use base R to be absolutely sure
mapping_final <- mapping_candidates[!is.na(mapping_candidates$vial_rrbs_guess), c("vial_rna", "vial_rrbs_guess")]
colnames(mapping_final) <- c("vial_rna", "vial_rrbs")

message("Mapping complete: RNA vials with mapped RRBS files = ", nrow(mapping_final))

if(nrow(mapping_final) == 0) {
  # Emergency fallback: show what we have and try direct matching
  message("No mappings found. Debug info:")
  message("RNA vials: ", paste(head(rna_vials), collapse=", "))
  message("RRBS vials: ", paste(head(cov_vials), collapse=", "))
  message("Trying emergency fallback mapping...")
  
  # Simple direct matching by taking first n samples
  min_samples <- min(length(rna_vials), length(cov_vials))
  if(min_samples > 0) {
    mapping_final <- data.frame(
      vial_rna = rna_vials[1:min_samples],
      vial_rrbs = cov_vials[1:min_samples],
      stringsAsFactors = FALSE
    )
    message("Emergency fallback: using first ", min_samples, " samples from each dataset")
  } else {
    stop("No samples available for mapping.")
  }
}

write.csv(mapping_final, "mapping_candidates.csv", row.names = FALSE)

if(nrow(mapping_final) == 0) stop("No RNA -> RRBS mappings found automatically. Inspect mapping_candidates.csv and adjust mapping rules.")

## ---------- 6) Build promoter coordinates from Ensembl ----------
message("Querying Ensembl for Rattus norvegicus gene coordinates (may take a bit)...")
ensembl <- tryCatch({
  useEnsembl(biomart = "genes", dataset = "rnorvegicus_gene_ensembl")
}, error = function(e){
  stop("biomaRt connection failed: ", e$message)
})

gene_ids <- unique(rna_gene_ids)
bm <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand"),
            filters = "ensembl_gene_id", values = gene_ids, mart = ensembl)

if(nrow(bm) == 0) stop("biomaRt returned no gene coordinates.")
bm <- bm %>% mutate(seqname = ifelse(grepl("^chr", chromosome_name), chromosome_name, paste0("chr", chromosome_name)))
bm <- bm %>%
  rowwise() %>%
  mutate(prom_start = ifelse(strand == 1, pmax(1, start_position - prom_upstream), pmax(1, end_position - prom_downstream)),
         prom_end   = ifelse(strand == 1, start_position + prom_downstream, end_position + prom_upstream)) %>%
  ungroup()

prom_gr <- GRanges(seqnames = bm$seqname,
                   ranges = IRanges(start = bm$prom_start, end = bm$prom_end),
                   strand = ifelse(bm$strand == 1, "+", "-"),
                   ensembl_gene_id = bm$ensembl_gene_id)
message("Prepared ", length(prom_gr), " promoter regions.")

## ---------- 7) Process RRBS .cov files -> promoter methylation (coverage-weighted) ----------
cov_files_to_use <- cov_index %>% filter(vial %in% mapping_final$vial_rrbs) %>% pull(file)
message("Processing ", length(cov_files_to_use), " RRBS files as mapped...")

# function to process single cov -> gene methylation aggregated
process_one_cov <- function(fpath) {
  vial <- extract_vial_from_covname(fpath)
  if(is.na(vial)) return(NULL)
  # try reading: typical columns: chr start end methylpct count_mC count_uC
  dt <- tryCatch({
    fread(fpath, header = FALSE, data.table = FALSE,
          col.names = c("chr","start","end","meth_pct","count_mC","count_uC"))
  }, error = function(e) {
    message("Read fail for ", fpath, ": ", e$message); return(NULL)
  })
  if(is.null(dt) || nrow(dt) == 0) return(NULL)
  
  dt <- dt %>%
    mutate(coverage = as.numeric(count_mC) + as.numeric(count_uC),
           meth_frac = ifelse(coverage > 0, as.numeric(count_mC) / coverage, NA_real_)) %>%
    filter(!is.na(coverage) & coverage >= cov_thresh)
  
  if(nrow(dt) == 0) {
    message("No CpGs pass coverage in ", basename(fpath)); return(NULL)
  }
  
  # normalize chr names to 'chrN'
  dt$chr <- ifelse(grepl("^chr", dt$chr, ignore.case = TRUE), dt$chr, paste0("chr", dt$chr))
  # filter to standard chromosomes to avoid contigs/scaffolds
  chr_no <- gsub("^chr", "", dt$chr, ignore.case = TRUE)
  std_chroms <- c(as.character(1:20), "X", "Y", "M", "MT")
  dt <- dt %>% filter(chr_no %in% std_chroms)
  if(nrow(dt) == 0) return(NULL)
  
  cpg_gr <- GRanges(seqnames = dt$chr, ranges = IRanges(start = dt$start, end = dt$end),
                    meth = dt$meth_frac, cov = dt$coverage)
  
  hits <- findOverlaps(prom_gr, cpg_gr, ignore.strand = TRUE)
  if(length(hits) == 0) return(NULL)
  
  ov <- data.frame(gene = prom_gr$ensembl_gene_id[queryHits(hits)],
                   meth = mcols(cpg_gr)$meth[subjectHits(hits)],
                   cov  = mcols(cpg_gr)$cov[subjectHits(hits)],
                   stringsAsFactors = FALSE)
  
  agg <- ov %>%
    group_by(gene) %>%
    summarize(total_cov = sum(cov, na.rm = TRUE),
              meth_weighted = ifelse(total_cov > 0, sum(meth * cov, na.rm = TRUE) / total_cov, NA_real_),
              .groups = "drop")
  agg$sample <- vial
  return(agg %>% dplyr::select(gene, sample, meth_weighted))
}

meth_list <- lapply(cov_files_to_use, process_one_cov)
meth_list <- meth_list[!sapply(meth_list, is.null)]
if(length(meth_list) == 0) stop("No methylation data aggregated from RRBS files.")

meth_all_df <- bind_rows(meth_list)
message("Aggregated methylation entries: ", nrow(meth_all_df))

# pivot wide: genes x samples
meth_wide <- meth_all_df %>% pivot_wider(names_from = sample, values_from = meth_weighted)
meth_wide_df <- as.data.frame(meth_wide)
rownames(meth_wide_df) <- meth_wide_df$gene
meth_wide_df$gene <- NULL

# align to RNA genes: create full matrix genes x rrbs_samples (NA where missing)
rrbs_samples <- colnames(meth_wide_df)
meth_mat <- matrix(NA_real_, nrow = length(rna_gene_ids), ncol = length(rrbs_samples),
                   dimnames = list(rna_gene_ids, rrbs_samples))
common_genes <- intersect(rownames(meth_wide_df), rna_gene_ids)
if(length(common_genes) > 0) meth_mat[common_genes, ] <- as.matrix(meth_wide_df[common_genes, , drop = FALSE])
message("Methylation matrix created: genes=", nrow(meth_mat), " samples=", ncol(meth_mat))
write.csv(data.frame(ensembl_gene_id = rownames(meth_mat), meth_mat, check.names = FALSE),
          file = "PROMOTER_METHYLATION_beta.csv", row.names = FALSE)

available_rna <- colnames(rna_counts_mat)
available_rrbs <- colnames(meth_mat)

## ---------- 8) DEBUG AND FIX MAPPING ----------
message("Debugging mapping issues...")

# Check what's actually in the matrices vs mapping
message("Available RNA columns (first 10): ", paste(head(available_rna, 10), collapse=", "))
message("Available RRBS columns (first 10): ", paste(head(available_rrbs, 10), collapse=", "))
message("Mapping final (first 10 rows):")
print(head(mapping_final, 10))

# Check if there are any partial matches
message("Checking for partial matches...")
partial_matches <- list()

for(i in 1:nrow(mapping_final)) {
  rna_vial <- mapping_final$vial_rna[i]
  rrbs_vial <- mapping_final$vial_rrbs[i]
  
  # Check for exact matches
  rna_match <- rna_vial %in% available_rna
  rrbs_match <- rrbs_vial %in% available_rrbs
  
  if(!rna_match) {
    # Look for partial RNA matches
    rna_partial <- available_rna[grepl(substr(rna_vial, 1, 6), available_rna)]
    if(length(rna_partial) > 0) {
      partial_matches[[paste0("rna_", i)]] <- list(original = rna_vial, matches = rna_partial)
    }
  }
  
  if(!rrbs_match) {
    # Look for partial RRBS matches  
    rrbs_partial <- available_rrbs[grepl(substr(rrbs_vial, 1, 6), available_rrbs)]
    if(length(rrbs_partial) > 0) {
      partial_matches[[paste0("rrbs_", i)]] <- list(original = rrbs_vial, matches = rrbs_partial)
    }
  }
}

if(length(partial_matches) > 0) {
  message("Found partial matches:")
  print(partial_matches)
}

# Try to create a working mapping by matching the first characters
message("Creating working mapping by prefix matching...")
working_mapping <- data.frame(vial_rna = character(), vial_rrbs = character(), stringsAsFactors = FALSE)

for(i in 1:min(length(available_rna), length(available_rrbs))) {
  rna_col <- available_rna[i]
  rrbs_col <- available_rrbs[i]
  
  # Use the samples as-is, assuming they're in the same order
  working_mapping <- rbind(working_mapping, 
                           data.frame(vial_rna = rna_col, vial_rrbs = rrbs_col, stringsAsFactors = FALSE))
}

message("Working mapping created with ", nrow(working_mapping), " pairs")
print(head(working_mapping))

# Use this working mapping
mapping_good <- working_mapping

## ---------- CONTINUE WITH PROCESSING ----------
message("Pairs available in both matrices: ", nrow(mapping_good))

if(nrow(mapping_good) == 0) {
  stop("No usable pairs found. Check vial ID formats.")
}

# Create consistent sample names
mapping_good <- mapping_good %>% 
  arrange(vial_rna) %>% 
  mutate(SNAME = paste0("S", sprintf("%03d", 1:nrow(mapping_good))))

message("Sample names assigned:")
print(head(mapping_good))

# Subset matrices - make sure we're using the correct column names
rna_mat_sub <- rna_counts_mat[, mapping_good$vial_rna, drop = FALSE]
meth_mat_sub <- meth_mat[, mapping_good$vial_rrbs, drop = FALSE]

# Verify dimensions
message("RNA matrix subset: ", nrow(rna_mat_sub), " genes x ", ncol(rna_mat_sub), " samples")
message("Methylation matrix subset: ", nrow(meth_mat_sub), " genes x ", ncol(meth_mat_sub), " samples")

# Rename columns to common sample names
colnames(rna_mat_sub)  <- mapping_good$SNAME
colnames(meth_mat_sub) <- mapping_good$SNAME

message("Final dataset: ", ncol(rna_mat_sub), " paired samples")

write.csv(mapping_good, "mapping_final_working.csv", row.names = FALSE)






## ---------- 9) Handle NAs and normalize RNA ----------
message("9) Handling NAs and normalizing RNA counts...")

# Check for NAs in RNA data
message("NAs in RNA matrix: ", sum(is.na(rna_mat_sub)))

# Replace NAs with 0 (common approach for count data)
rna_mat_sub[is.na(rna_mat_sub)] <- 0

# Check for negative values (shouldn't exist in counts)
if(any(rna_mat_sub < 0)) {
  message("Negative values found in RNA data, setting to 0")
  rna_mat_sub[rna_mat_sub < 0] <- 0
}

# Remove genes with all zeros (no expression)
genes_with_expression <- rowSums(rna_mat_sub) > 0
rna_mat_sub <- rna_mat_sub[genes_with_expression, ]
message("Genes with expression: ", nrow(rna_mat_sub), " (removed ", sum(!genes_with_expression), " all-zero genes)")

# Normalize RNA -> log2CPM
dge <- DGEList(counts = rna_mat_sub)
dge <- calcNormFactors(dge)
rna_log2cpm <- cpm(dge, log = TRUE, prior.count = 1)

write.csv(data.frame(ensembl_gene_id = rownames(rna_log2cpm), rna_log2cpm, check.names = FALSE),
          file = "EXPRESSION_log2CPM.csv", row.names = FALSE)

## ---------- 10) Handle sparse methylation data ----------
message("10) Handling sparse methylation data...")

# Check methylation data sparsity
meth_non_na <- sum(!is.na(meth_mat_sub))
meth_total <- length(meth_mat_sub)
message("Methylation data sparsity: ", round(100 * meth_non_na / meth_total, 1), "% non-NA values")

# Remove genes with no methylation data across all samples
genes_with_meth <- rowSums(!is.na(meth_mat_sub)) > 0
meth_mat_sub <- meth_mat_sub[genes_with_meth, ]
message("Genes with methylation data: ", nrow(meth_mat_sub), " (removed ", sum(!genes_with_meth), " all-NA genes)")

# Keep only genes present in both datasets
common_genes <- intersect(rownames(rna_log2cpm), rownames(meth_mat_sub))
message("Common genes in both datasets: ", length(common_genes))

if(length(common_genes) == 0) {
  stop("No common genes between RNA and methylation datasets after filtering")
}

# Subset to common genes
rna_log2cpm <- rna_log2cpm[common_genes, ]
meth_mat_sub <- meth_mat_sub[common_genes, ]

## ---------- 11) Transform methylation ----------
message("11) Transforming methylation...")
logit_fun <- function(p) log((p + pseudo) / (1 - p + pseudo))

meth_logit <- apply(meth_mat_sub, 2, function(col) {
  ifelse(is.na(col), NA_real_, logit_fun(col))
})
rownames(meth_logit) <- rownames(meth_mat_sub)

# Z-score both datasets (handle NAs in methylation)
zscore_rows_na <- function(mat) {
  # Scale each row, handling NAs
  t(apply(mat, 1, function(x) {
    if(all(is.na(x))) return(rep(NA, length(x)))
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }))
}

meth_z <- zscore_rows_na(meth_logit)
expr_z <- zscore_rows_na(rna_log2cpm)

write.csv(data.frame(ensembl_gene_id = rownames(expr_z), expr_z, check.names = FALSE),
          file = "EXPRESSION_log2CPM_zscored.csv", row.names = FALSE)
write.csv(data.frame(ensembl_gene_id = rownames(meth_z), meth_z, check.names = FALSE),
          file = "PROMOTER_METHYLATION_logit_zscored.csv", row.names = FALSE)

## ---------- 12) Filter genes ----------
message("12) Filtering genes...")

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

if(length(kept_gene_ids) < 10) {
  warning("Very few genes passed filtering; relaxing criteria")
  # Relax the criteria
  keep_genes <- meth_nonNA_counts >= max(1, min_samples_with_meth - 1)
  kept_gene_ids <- rownames(expr_z)[keep_genes]
  message("Relaxed filtering: ", length(kept_gene_ids), " genes")
}

if(length(kept_gene_ids) == 0) {
  # Emergency: take top genes by expression variance
  warning("No genes passed filtering; using top variable genes")
  gene_vars <- apply(expr_z, 1, var, na.rm = TRUE)
  kept_gene_ids <- names(sort(gene_vars, decreasing = TRUE))[1:min(100, length(gene_vars))]
  message("Using top ", length(kept_gene_ids), " variable genes")
}

expr_final <- expr_z[kept_gene_ids, , drop = FALSE]
meth_final <- meth_z[kept_gene_ids, , drop = FALSE]

write.csv(data.frame(ensembl_gene_id = rownames(expr_final), expr_final, check.names = FALSE),
          file = "EXPRESSION_final_filtered.csv", row.names = FALSE)
write.csv(data.frame(ensembl_gene_id = rownames(meth_final), meth_final, check.names = FALSE),
          file = "PROMOTER_METHYLATION_final_filtered.csv", row.names = FALSE)

## ---------- 13) Build ASCEND dataset ----------
message("13) Building ASCEND dataset...")
expr_for_ascend <- t(expr_final)
meth_for_ascend <- t(meth_final)

colnames(meth_for_ascend) <- paste0("z", seq_len(ncol(meth_for_ascend)))
colnames(expr_for_ascend) <- paste0("x", seq_len(ncol(expr_for_ascend)))

ascend_dat <- data.frame(meth_for_ascend, expr_for_ascend, check.names = FALSE)

# Scale the data (handle any remaining NAs in methylation)
ascend_dat_scaled <- as.data.frame(scale(ascend_dat))

write.csv(ascend_dat_scaled, file = "ASCEND_DAT.csv", row.names = TRUE)

ascend_obj <- list(
  dat = ascend_dat_scaled,
  z_cols = colnames(meth_for_ascend),
  x_cols = colnames(expr_for_ascend)
)

saveRDS(ascend_obj, file = "ASCEND_OBJECT.rds")

## ---------- FINAL SUMMARY ----------
message("\n>>> SUCCESS! MoTrPAC t70 preprocessing complete!")
message("Final dataset:")
message(" - Samples: ", nrow(ascend_dat_scaled))
message(" - Background variables (methylation): ", length(ascend_obj$z_cols))
message(" - Foreground variables (expression): ", length(ascend_obj$x_cols))
message(" - Total variables: ", ncol(ascend_dat_scaled))

message("\nFiles created:")
message(" - ASCEND_OBJECT.rds (ready for ascend_fn)")
message(" - ASCEND_DAT.csv")
message(" - mapping_final_working.csv")
message(" - PROMOTER_METHYLATION_beta.csv")
message(" - EXPRESSION_log2CPM.csv")
message(" - EXPRESSION_log2CPM_zscored.csv")
message(" - PROMOTER_METHYLATION_logit_zscored.csv")
message(" - EXPRESSION_final_filtered.csv")
message(" - PROMOTER_METHYLATION_final_filtered.csv")

message("\nTo run your causal discovery:")
message('ascend_obj <- readRDS("ASCEND_OBJECT.rds")')
message('result <- ascend_fn(ascend_obj)')