## ======================================================================
## tests/make_fake_dgrp.R
## Writes small fake DGRP files in the real files' formats (genotype table
## with snp_id and line_* columns coded 0/2/-, GEO expression table with
## "Gene ID" and *_M_mean columns), with planted SNP -> gene -> gene effects,
## so the DGRP control pipeline can be tested without the real data.
##
## Reviewer comments: 7 (tests the DGRP control scripts)
## How to run: Rscript tests/make_fake_dgrp.R <out_dir> [n_snps] [n_genes]
## ======================================================================
library(data.table)
a <- commandArgs(trailingOnly = TRUE)
out <- a[1]; n_snp <- if (length(a) > 1) as.integer(a[2]) else 3000L
n_gene <- if (length(a) > 2) as.integer(a[3]) else 2000L
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(7)
lines <- sort(sample(20:900, 205)); n <- length(lines)
## genotypes: 0/2 with allele frequencies 0.05-0.5, 5% of SNPs with a missing call
af <- runif(n_snp, 0.05, 0.5)
G <- t(sapply(af, function(p) 2 * rbinom(n, 1, p)))
geno <- as.data.table(matrix(as.character(G), n_snp))
setnames(geno, paste0("line_", lines))
miss <- sample(n_snp, round(0.05 * n_snp))
for (i in miss) set(geno, i, sample(n, 1), "-")
geno <- cbind(data.table(snp_id = sprintf("2R_%d_SNP", seq_len(n_snp))), geno)
fwrite(geno, file.path(out, "dgrp2.tgeno.txt"), sep = "\t")
## expression: 300 variable genes driven by common SNPs; 5 hub genes drive 30 others
X <- matrix(rnorm(n_gene * n, 7, 0.3), n_gene)
common <- which(af > 0.25)
var_genes <- 1:300
for (g in var_genes) {
  s <- sample(common, 2)
  X[g, ] <- 7 + 0.6 * scale(G[s[1], ])[, 1] + 0.4 * scale(G[s[2], ])[, 1] + rnorm(n, 0, 0.6)
}
hubs <- 1:5
for (h in hubs) for (g in sample(6:300, 30)) X[g, ] <- X[g, ] + 0.5 * (X[h, ] - 7)
X[sample(length(X), 50)] <- NA      # a few genes with missing values
expr <- data.table(`Gene ID` = sprintf("FBgn%07d", seq_len(n_gene)))
expr <- cbind(expr, as.data.table(X))
setnames(expr, c("Gene ID", paste0(lines, "_M_mean")))
fwrite(expr, file.path(out, "GSE117850_DGRP_GEO_Table_4_Gene_Male_Line_Means.txt.gz"), sep = "\t")
writeLines(sprintf("FBgn%07d", c(1:5, sample(6:2000, 40))), file.path(out, "immune_genes.txt"))
cat("fake DGRP data in", out, "\n")
