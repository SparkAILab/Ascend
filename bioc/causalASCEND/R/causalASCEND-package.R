#' causalASCEND: scalable two-tier causal discovery
#'
#' Causal ASCEND infers ancestral relations among foreground variables
#' (for example gene expression) when a background layer (for example
#' genotypes) is known to precede them. Start with [ascend()] (matrices or
#' `SummarizedExperiment` objects) or [ascendMultiAssay()]
#' (`MultiAssayExperiment`), and read the result with [ancestralEdges()].
#' The vignette walks through a complete analysis:
#' `vignette("causalASCEND")`.
#'
#' @keywords internal
#' @importFrom methods new is slot setClass setGeneric setMethod show
#'   setValidity
#' @importFrom utils head
#' @importFrom SummarizedExperiment assay
"_PACKAGE"
