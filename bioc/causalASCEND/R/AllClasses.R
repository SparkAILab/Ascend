#' The result of an ASCEND run
#'
#' `AscendResult` holds the ancestral relations that [ascend()] inferred
#' among the foreground variables, together with the evidence behind each
#' one and an account of the work done. Use the accessors listed below
#' rather than the slots.
#'
#' @slot ancestral Numeric matrix, foreground by foreground. Entry `[a, b]`
#'   is 1 if `a` is an ancestor of `b`, 0.5 if `a` is not a descendant of
#'   `b`, 0 if neither is an ancestor of the other, and `NA` if the pair is
#'   unresolved. Rows and columns follow a topological order of the claims.
#' @slot pValues Numeric matrix: p-value of the last pairwise
#'   (rule R3) independence test of each pair.
#' @slot rule Character matrix: the step that committed each pair (`"R1"`,
#'   `"R2"`, `"R3"`, `"transitivity"`, `"symmetry"`, `"R3_final"`), or
#'   `NA` if unresolved.
#' @slot ruleSweep Integer matrix: the sweep in which the pair was
#'   committed.
#' @slot ruleVotes Integer matrix: for R1 and R2 commits, the number of
#'   agreeing witnesses.
#' @slot stats List: numbers of conditional independence tests by type,
#'   time spent in them, conditioning-set sizes and number of sweeps.
#' @slot params List: the arguments of the call and the data dimensions.
#'
#' @seealso [ascend()], [ancestralMatrix()], [ancestralEdges()],
#'   [pairPValues()], [ascendStats()]
#' @return An object of class `AscendResult`.
#' @examples
#' set.seed(1)
#' sim <- simulateTwoTier(n = 500, nBackground = 10, nForeground = 6)
#' res <- ascend(sim$background, sim$foreground)
#' res
#' @exportClass AscendResult
setClass("AscendResult",
    slots = c(
        ancestral = "matrix",
        pValues = "matrix",
        rule = "matrix",
        ruleSweep = "matrix",
        ruleVotes = "matrix",
        stats = "list",
        params = "list"))

setValidity("AscendResult", function(object) {
    M <- object@ancestral
    msg <- character(0)
    if (nrow(M) != ncol(M))
        msg <- c(msg, "'ancestral' must be a square matrix")
    if (!identical(rownames(M), colnames(M)))
        msg <- c(msg, "'ancestral' must have identical row and column names")
    ok <- is.na(M) | M %in% c(0, 0.5, 1)
    if (!all(ok))
        msg <- c(msg, "'ancestral' entries must be 1, 0.5, 0 or NA")
    for (s in c("pValues", "rule", "ruleSweep", "ruleVotes")) {
        if (!identical(dim(methods::slot(object, s)), dim(M)))
            msg <- c(msg, sprintf(
                "'%s' must have the dimensions of 'ancestral'", s))
    }
    if (length(msg)) msg else TRUE
})
