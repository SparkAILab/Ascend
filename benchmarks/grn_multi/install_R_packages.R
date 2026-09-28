## ======================================================================
## benchmarks/grn_multi/install_R_packages.R
## Installs every R package the benchmark needs (CRAN + Bioconductor).
## How to run (once, on a login node or in an interactive job):
##   Rscript benchmarks/grn_multi/install_R_packages.R
## Installs into the first writable library in .libPaths(); set R_LIBS_USER
## first if you want a specific location.
## ======================================================================
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4)
cran <- c("data.table", "jsonlite", "dplyr", "tidyverse", "foreach", "doMC", "glmnet",
          "lightgbm", "matrixStats", "R.utils", "PRROC", "ppcor", "bnlearn", "WGCNA",
          "ggplot2", "igraph", "BiocManager")
bioc <- c("graph", "RBGL", "minet", "GENIE3", "impute", "preprocessCore", "GO.db")
need <- setdiff(cran, rownames(installed.packages()))
if (length(need)) install.packages(need)
need <- setdiff(bioc, rownames(installed.packages()))
if (length(need)) BiocManager::install(need, update = FALSE, ask = FALSE)
for (p in c("pcalg", "tpc")) if (!p %in% rownames(installed.packages())) install.packages(p)
all <- c(cran, bioc, "pcalg", "tpc")
ok <- vapply(all, requireNamespace, logical(1), quietly = TRUE)
print(ok)
if (!all(ok)) stop("Not installed: ", paste(all[!ok], collapse = ", "))
cat("All R packages installed.\n")
