#' @include AllClasses.R
NULL

#' Accessors for AscendResult objects
#'
#' Extract the inferred ancestral relations and the evidence behind them
#' from the result of [ascend()].
#'
#' @param x An [AscendResult-class] object.
#' @param include Which relations `ancestralEdges()` returns: any of
#'   `"ancestor"`, `"non-descendant"`, `"unrelated"` and `"unresolved"`.
#' @param object An [AscendResult-class] object.
#'
#' @return
#' * `ancestralMatrix()`: the numeric ancestrality matrix (1 ancestor,
#'   0.5 non-descendant, 0 unrelated, `NA` unresolved), row variable
#'   relative to column variable.
#' * `pairPValues()`: matrix of the p-value of the last pairwise
#'   independence test of each pair, a continuous score for ranking pairs.
#' * `ascendStats()`: list with the numbers of conditional independence
#'   tests by type (`nPairTests`, `nWitnessTests`, `nBlanketTests`), the
#'   time spent in them, conditioning-set sizes and the number of sweeps.
#' * `ancestralEdges()`: data frame with one row per relation, with
#'   columns `from`, `to`, `relation`, `rule`, `sweep`, `votes` and
#'   `pValue`. Ancestor and non-descendant relations are listed once, from
#'   the earlier to the later variable; unrelated and unresolved pairs are
#'   listed once with `from` before `to` in the matrix order.
#'
#' @examples
#' set.seed(1)
#' sim <- simulateTwoTier(n = 500, nBackground = 10, nForeground = 6)
#' res <- ascend(sim$background, sim$foreground)
#' ancestralMatrix(res)
#' ancestralEdges(res)
#' ascendStats(res)$nPairTests
#' @name AscendResult-accessors
NULL

#' @rdname AscendResult-accessors
#' @export
setGeneric("ancestralMatrix", function(x) standardGeneric("ancestralMatrix"))

#' @rdname AscendResult-accessors
#' @export
setMethod("ancestralMatrix", "AscendResult", function(x) x@ancestral)

#' @rdname AscendResult-accessors
#' @export
setGeneric("pairPValues", function(x) standardGeneric("pairPValues"))

#' @rdname AscendResult-accessors
#' @export
setMethod("pairPValues", "AscendResult", function(x) x@pValues)

#' @rdname AscendResult-accessors
#' @export
setGeneric("ascendStats", function(x) standardGeneric("ascendStats"))

#' @rdname AscendResult-accessors
#' @export
setMethod("ascendStats", "AscendResult", function(x) x@stats)

#' @rdname AscendResult-accessors
#' @export
setGeneric("ancestralEdges",
    function(x, include = c("ancestor", "non-descendant"))
        standardGeneric("ancestralEdges"))

#' @rdname AscendResult-accessors
#' @export
setMethod("ancestralEdges", "AscendResult",
    function(x, include = c("ancestor", "non-descendant")) {
        include <- match.arg(include,
            c("ancestor", "non-descendant", "unrelated", "unresolved"),
            several.ok = TRUE)
        M <- x@ancestral
        nm <- rownames(M)
        rows <- list()
        for (a in seq_len(nrow(M))) for (b in seq_len(ncol(M))) {
            if (a == b) next
            v <- M[a, b]
            w <- M[b, a]
            rel <- if (is.na(v)) {
                if (a < b && is.na(w)) "unresolved" else NA
            } else if (v == 1) "ancestor"
            else if (v == 0.5) "non-descendant"
            else if (a < b && !is.na(w) && w == 0) "unrelated"
            else NA
            if (is.na(rel) || !rel %in% include) next
            rows[[length(rows) + 1]] <- data.frame(from = nm[a], to = nm[b],
                relation = rel, rule = x@rule[a, b],
                sweep = x@ruleSweep[a, b], votes = x@ruleVotes[a, b],
                pValue = x@pValues[a, b], stringsAsFactors = FALSE)
        }
        if (!length(rows))
            return(data.frame(from = character(0), to = character(0),
                relation = character(0), rule = character(0),
                sweep = integer(0), votes = integer(0), pValue = numeric(0)))
        out <- do.call(rbind, rows)
        rownames(out) <- NULL
        out
    })

#' @rdname AscendResult-accessors
#' @export
setMethod("show", "AscendResult", function(object) {
    M <- object@ancestral
    up <- M[upper.tri(M)]
    lo <- t(M)[upper.tri(M)]
    nPairs <- length(up)
    nAnc <- sum(up == 1 | lo == 1, na.rm = TRUE)
    nNd <- sum((up == 0.5 | lo == 0.5) & !(up == 1 | lo == 1), na.rm = TRUE)
    nUnrel <- sum(up == 0 & lo == 0, na.rm = TRUE)
    nNA <- sum(is.na(up) & is.na(lo))
    p <- object@params
    cat("AscendResult:", p$nForeground, "foreground and", p$nBackground,
        "background variables,", p$nSamples, "samples\n")
    cat(sprintf("  %d pairs: %d ancestor, %d non-descendant, %d unrelated,",
        nPairs, nAnc, nNd, nUnrel), nNA, "unresolved\n")
    cat(sprintf("  %d CI tests in %d sweep(s)\n", object@stats$nTests,
        object@stats$nSweeps))
    cat("  accessors: ancestralMatrix(), ancestralEdges(), pairPValues(),",
        "ascendStats()\n")
    invisible(object)
})
