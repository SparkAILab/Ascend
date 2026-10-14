# causalASCEND

Causal ASCEND (Ancestral Scalable Causal discovEry via iNherited Descent)
infers ancestral (cause-before-effect) relations among downstream
molecular variables, such as gene expression, when a layer of upstream
variables, such as genotypes, is known to precede them.

```r
BiocManager::install("causalASCEND")   # once accepted by Bioconductor
library(causalASCEND)
set.seed(1)
sim <- simulateTwoTier(n = 2000, nBackground = 20, nForeground = 15,
    r2 = 0.7, sparsity = 0.9, pCross = 0.2, effect = 0.9)
res <- ascend(sim$background, sim$foreground)
ancestralEdges(res)
```

Input: two matrices (features in rows, samples in columns), two
`SummarizedExperiment` objects, or a `MultiAssayExperiment`
(`ascendMultiAssay()`). See `vignette("causalASCEND")`.

Reference: Asiedu S, Watson D. Causal ASCEND: scalable two-tier causal
discovery on high-dimensional multi-omics data.
