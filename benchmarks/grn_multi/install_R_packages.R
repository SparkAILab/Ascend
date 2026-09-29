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
          "lightgbm", "matrixStats", "R.utils", "PRROC", "ppcor", "bnlearn",
          "ggplot2", "igraph", "BiocManager")
bioc <- c("graph", "RBGL", "minet", "GENIE3", "impute", "preprocessCore", "GO.db", "AnnotationDbi")
have <- function() rownames(installed.packages())
# one package at a time, so one failure does not stop the others and its
# compiler error is printed right above the next package
inst <- function(pkgs, fun) for (p in setdiff(pkgs, have())) {
  cat("\n==== installing", p, "\n")
  try(fun(p))
}
inst("BiocManager", install.packages)
inst(cran, install.packages)
inst(bioc, function(p) BiocManager::install(p, update = FALSE, ask = FALSE))
# WGCNA depends on the Bioconductor packages above, so it goes last
inst(c("WGCNA", "pcalg", "tpc"), install.packages)
cran <- c(cran, "WGCNA")
all <- c(cran, bioc, "pcalg", "tpc")
ok <- vapply(all, requireNamespace, logical(1), quietly = TRUE)
print(ok)
if (!all(ok)) stop("Not installed: ", paste(all[!ok], collapse = ", "),
                   ". Scroll up to its ==== line for the compiler error.")
cat("All R packages installed.\n")
