# causalASCEND 0.99.0

* Initial Bioconductor submission.
* `ascend()` for matrices, data frames and `SummarizedExperiment` objects;
  `ascendMultiAssay()` for `MultiAssayExperiment` objects.
* `AscendResult` class with `ancestralMatrix()`, `ancestralEdges()`,
  `pairPValues()` and `ascendStats()`.
* `simulateTwoTier()`, `trueAncestral()` and `evaluateAncestral()` for
  simulation and benchmarking.
* Results are identical to the implementation used in the paper.
