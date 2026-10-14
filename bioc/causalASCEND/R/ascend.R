#' @include AllClasses.R
NULL

#' Infer ancestral relations with Causal ASCEND
#'
#' `ascend()` infers which foreground variables (for example genes) are
#' causal ancestors of which, using a background layer (for example
#' genotypes) that is known to precede the foreground. It tests
#' conditional independence with a Fisher-z partial-correlation test,
#' conditioning each test only on the nearest ancestors of the two
#' variables, which it learns and refines over successive sweeps.
#'
#' The method assumes that no foreground variable causes a background
#' variable, and that the data are Markov and faithful to a directed
#' acyclic graph. Not every background variable needs to be measured, but
#' missing background variables act as hidden confounders and lower
#' accuracy. Missing values are handled with pairwise-complete
#' correlations.
#'
#' Each pair of foreground variables is resolved by one of three rules,
#' with background or already ordered foreground variables as witnesses:
#' R1 (deactivation) gives an ancestor relation, R2 (activation) a
#' non-descendant relation and R3 (conditional independence) an
#' "unrelated" call. Relations are then closed transitively, the nearest
#' ancestors are re-learned and the sweep repeats until nothing changes.
#'
#' @param background Background (upstream) data: a numeric matrix or
#'   data frame with features in rows and samples in columns, or a
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#' @param foreground Foreground (downstream) data, in the same form as
#'   `background` and with the same samples.
#' @param alpha Significance level of the orientation rules.
#' @param alphaMB Significance level of the Markov-blanket (nearest
#'   ancestor) searches.
#' @param fdr Logical: correct the pairwise tests of each sweep with
#'   the Benjamini-Hochberg procedure.
#' @param minVotes Number of agreeing witnesses required before a
#'   direction is committed. Values above 1 trade coverage for precision.
#' @param prescreen Benjamini-Hochberg threshold of the marginal
#'   prescreen that limits each foreground variable's background
#'   candidates. Set to 1 to keep every background variable.
#' @param maxIter Maximum number of sweeps.
#' @param verbose Logical: report progress with [message()].
#' @param backgroundAssay,foregroundAssay Assay name or index used when
#'   `background` or `foreground` is a `SummarizedExperiment`.
#'
#' @return An [AscendResult-class] object. See [ancestralMatrix()],
#'   [ancestralEdges()], [pairPValues()] and [ascendStats()].
#'
#' @references Asiedu S, Watson D. Causal ASCEND: scalable two-tier causal
#'   discovery on high-dimensional multi-omics data.
#'
#' @seealso [ascendMultiAssay()] for a `MultiAssayExperiment`,
#'   [simulateTwoTier()] for test data.
#' @examples
#' set.seed(123)
#' sim <- simulateTwoTier(n = 1000, nBackground = 20, nForeground = 8,
#'     r2 = 0.7, sparsity = 0.7, pCross = 0.1, effect = 1.5)
#' res <- ascend(sim$background, sim$foreground)
#' res
#' head(ancestralEdges(res))
#' evaluateAncestral(res, trueAncestral(sim))
#' @export
ascend <- function(background, foreground, alpha = 0.05, alphaMB = 0.05,
    fdr = TRUE, minVotes = 1, prescreen = 0.30, maxIter = 10,
    verbose = FALSE, backgroundAssay = 1, foregroundAssay = 1) {
    .checkScalar(alpha, "alpha", 0, 1)
    .checkScalar(alphaMB, "alphaMB", 0, 1)
    .checkScalar(prescreen, "prescreen", 0, 1)
    .checkScalar(minVotes, "minVotes", 1, Inf)
    .checkScalar(maxIter, "maxIter", 1, Inf)
    if (!is.logical(fdr) || length(fdr) != 1 || is.na(fdr))
        stop("'fdr' must be TRUE or FALSE")
    if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose))
        stop("'verbose' must be TRUE or FALSE")

    Z <- .asFeatureMatrix(background, backgroundAssay, "background", "z")
    X <- .asFeatureMatrix(foreground, foregroundAssay, "foreground", "x")
    X <- .alignSamples(Z, X)
    if (nrow(X) < 2) stop("'foreground' needs at least two features")
    if (nrow(Z) < 1) stop("'background' needs at least one feature")
    dup <- intersect(rownames(Z), rownames(X))
    if (length(dup))
        stop("feature names must differ between background and foreground; ",
            "shared: ", paste(utils::head(dup, 5), collapse = ", "))
    n <- ncol(X)
    if (n < 10) stop("at least 10 samples are needed")

    dat <- t(rbind(Z, X))
    dat[!is.finite(dat)] <- NA
    sds <- apply(dat, 2, stats::sd, na.rm = TRUE)
    const <- colnames(dat)[is.na(sds) | sds == 0]
    if (length(const))
        warning(length(const), " constant feature(s) treated as independent ",
            "of all others: ", paste(utils::head(const, 5), collapse = ", "))
    keep <- !colnames(dat) %in% const
    C <- diag(ncol(dat))
    dimnames(C) <- list(colnames(dat), colnames(dat))
    if (any(keep))
        C[keep, keep] <- stats::cor(dat[, keep, drop = FALSE],
            use = "pairwise.complete.obs")
    C[is.na(C)] <- 0
    diag(C) <- 1

    out <- .ascendEngine(C, n, rownames(Z), rownames(X), maxIter = maxIter,
        alpha = alpha, alphaMB = alphaMB, fdr = fdr, minVotes = minVotes,
        prescreen = prescreen, verbose = verbose)
    methods::new("AscendResult",
        ancestral = out$M, pValues = out$pairP, rule = out$rule,
        ruleSweep = out$ruleSweep, ruleVotes = out$ruleVotes,
        stats = out$stats,
        params = list(alpha = alpha, alphaMB = alphaMB, fdr = fdr,
            minVotes = minVotes, prescreen = prescreen, maxIter = maxIter,
            nSamples = n, nBackground = nrow(Z), nForeground = nrow(X)))
}

#' Run ASCEND on two experiments of a MultiAssayExperiment
#'
#' Extracts the background and foreground experiments, keeps the samples
#' (primary units) present in both and passes them to [ascend()].
#'
#' @param mae A
#'   \code{\link[MultiAssayExperiment]{MultiAssayExperiment}} object.
#' @param background,foreground Names (or indices) of the background and
#'   foreground experiments in `mae`.
#' @param backgroundAssay,foregroundAssay Assay used when an experiment
#'   is a `SummarizedExperiment`.
#' @param ... Further arguments passed to [ascend()].
#' @return An [AscendResult-class] object.
#' @examples
#' if (requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
#'     set.seed(2)
#'     sim <- simulateTwoTier(n = 300, nBackground = 10, nForeground = 5)
#'     ids <- paste0("s", seq_len(300))
#'     colnames(sim$background) <- colnames(sim$foreground) <- ids
#'     mae <- MultiAssayExperiment::MultiAssayExperiment(
#'         experiments = list(snp = sim$background, rna = sim$foreground))
#'     ascendMultiAssay(mae, background = "snp", foreground = "rna")
#' }
#' @export
ascendMultiAssay <- function(mae, background, foreground,
    backgroundAssay = 1, foregroundAssay = 1, ...) {
    if (!requireNamespace("MultiAssayExperiment", quietly = TRUE))
        stop("ascendMultiAssay() needs the MultiAssayExperiment package")
    if (!methods::is(mae, "MultiAssayExperiment"))
        stop("'mae' must be a MultiAssayExperiment")
    exps <- MultiAssayExperiment::experiments(mae)
    smap <- as.data.frame(MultiAssayExperiment::sampleMap(mae))
    pick <- function(which, arg) {
        nm <- if (is.numeric(which)) names(exps)[which] else which
        if (length(nm) != 1 || is.na(nm) || !nm %in% names(exps))
            stop("'", arg, "' must name one experiment of 'mae'")
        nm
    }
    bName <- pick(background, "background")
    fName <- pick(foreground, "foreground")
    toPrimary <- function(nm) {
        m <- smap[smap$assay == nm, , drop = FALSE]
        stats::setNames(as.character(m$primary), as.character(m$colname))
    }
    bMap <- toPrimary(bName)
    fMap <- toPrimary(fName)
    bExp <- exps[[bName]]
    fExp <- exps[[fName]]
    bPrim <- bMap[colnames(bExp)]
    fPrim <- fMap[colnames(fExp)]
    common <- intersect(bPrim[!is.na(bPrim)], fPrim[!is.na(fPrim)])
    common <- common[!common %in% c(bPrim[duplicated(bPrim)],
        fPrim[duplicated(fPrim)])]
    if (length(common) < 10)
        stop("fewer than 10 samples are shared by the two experiments")
    bExp <- bExp[, match(common, bPrim)]
    fExp <- fExp[, match(common, fPrim)]
    colnames(bExp) <- colnames(fExp) <- common
    ascend(bExp, fExp, backgroundAssay = backgroundAssay,
        foregroundAssay = foregroundAssay, ...)
}

## ---- input helpers ----------------------------------------------------

.checkScalar <- function(x, name, lo, hi) {
    if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < lo || x > hi)
        stop(sprintf("'%s' must be a single number in [%s, %s]", name, lo,
            hi))
    invisible(TRUE)
}

## Return a numeric matrix with features in rows and samples in columns.
.asFeatureMatrix <- function(x, assay, arg, prefix) {
    if (methods::is(x, "SummarizedExperiment"))
        x <- SummarizedExperiment::assay(x, assay)
    if (is.data.frame(x)) x <- as.matrix(x)
    if (methods::is(x, "Matrix") || methods::is(x, "DelayedArray"))
        x <- as.matrix(x)
    if (!is.matrix(x) || !is.numeric(x))
        stop("'", arg, "' must be a numeric matrix, data frame or ",
            "SummarizedExperiment with features in rows and samples in ",
            "columns")
    storage.mode(x) <- "double"
    if (is.null(rownames(x)))
        rownames(x) <- paste0(prefix, seq_len(nrow(x)))
    if (anyDuplicated(rownames(x)))
        stop("feature (row) names of '", arg, "' must be unique")
    x
}

## Put the foreground samples in the order of the background samples.
.alignSamples <- function(Z, X) {
    if (!is.null(colnames(Z)) && !is.null(colnames(X))) {
        if (!setequal(colnames(Z), colnames(X)) || anyDuplicated(colnames(Z)))
            stop("'background' and 'foreground' must have the same, unique ",
                "sample (column) names")
        return(X[, colnames(Z), drop = FALSE])
    }
    if (ncol(Z) != ncol(X))
        stop("'background' and 'foreground' must have the same number of ",
            "samples (columns)")
    X
}
