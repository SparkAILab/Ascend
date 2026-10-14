#' @include AllClasses.R
NULL

#' True ancestral relations of a simulated graph
#'
#' @param x The output of [simulateTwoTier()], or a 0/1 adjacency matrix
#'   with `x[parent, child] = 1`.
#' @return A 0/1 matrix with entry `[a, b] = 1` if `a` is an ancestor of
#'   `b` (diagonal `NA`), in a topological order.
#' @examples
#' set.seed(3)
#' sim <- simulateTwoTier(n = 100, nBackground = 5, nForeground = 5,
#'     sparsity = 0.5)
#' trueAncestral(sim)
#' @export
trueAncestral <- function(x) {
    A <- if (is.list(x) && !is.null(x$adjacency)) x$adjacency else x
    if (!is.matrix(A) || nrow(A) != ncol(A))
        stop("'x' must be simulateTwoTier() output or a square adjacency ",
            "matrix")
    A[is.na(A)] <- 0
    A <- (A != 0) * 1
    diag(A) <- 0
    d <- nrow(A)
    reach <- A
    P <- A
    for (k in seq_len(max(d - 1, 0))) {
        P <- ((P %*% A) > 0) * 1
        reach <- ((reach + P) > 0) * 1
    }
    dimnames(reach) <- dimnames(A)
    ord <- .topoOrder(reach)
    if (!is.null(ord)) reach <- reach[ord, ord, drop = FALSE]
    diag(reach) <- NA
    reach
}

#' Score inferred ancestral relations against the truth
#'
#' Each ordered pair `(a, b)` is a claim that `a` precedes `b` (entry 1
#' or 0.5) or not (entry 0). A reversed claim therefore counts as both a
#' false positive and a false negative. Unresolved pairs are counted
#' separately and are not scored.
#'
#' @param estimate An [AscendResult-class] object or an ancestrality
#'   matrix.
#' @param truth The output of [trueAncestral()].
#' @return A list with `tp`, `fp`, `fn`, `tn`, `unresolved`, `precision`,
#'   `recall`, `f1`, `orientationAccuracy` (fraction of claims on truly
#'   related pairs that point the right way) and `coverage` (fraction of
#'   pairs with any verdict).
#' @examples
#' set.seed(123)
#' sim <- simulateTwoTier(n = 1000, nBackground = 20, nForeground = 8,
#'     r2 = 0.7, sparsity = 0.7, pCross = 0.1, effect = 1.5)
#' res <- ascend(sim$background, sim$foreground)
#' evaluateAncestral(res, trueAncestral(sim))
#' @export
evaluateAncestral <- function(estimate, truth) {
    est <- if (methods::is(estimate, "AscendResult"))
        ancestralMatrix(estimate) else estimate
    nms <- rownames(truth)
    if (is.null(nms) || !setequal(nms, rownames(est)))
        stop("'estimate' and 'truth' must have the same variable names")
    est <- est[nms, nms, drop = FALSE]
    d <- nrow(truth)
    isPos <- function(v) !is.na(v) & v %in% c(0.5, 1)
    tp <- fp <- fn <- tn <- unresolved <- 0
    for (a in seq_len(d)) for (b in seq_len(d)) {
        if (a == b) next
        t <- !is.na(truth[a, b]) && truth[a, b] == 1
        e <- est[a, b]
        if (is.na(e)) {
            unresolved <- unresolved + 1
            next
        }
        p <- isPos(e)
        if (t && p) tp <- tp + 1
        else if (!t && p) fp <- fp + 1
        else if (t && !p) fn <- fn + 1
        else tn <- tn + 1
    }
    dirOk <- dirTot <- resolved <- 0
    pairs <- if (d > 1) which(upper.tri(truth), arr.ind = TRUE) else
        matrix(0L, 0, 2)
    for (r in seq_len(nrow(pairs))) {
        a <- pairs[r, 1]
        b <- pairs[r, 2]
        if (!is.na(est[a, b]) || !is.na(est[b, a])) resolved <- resolved + 1
        tab <- !is.na(truth[a, b]) && truth[a, b] == 1
        tba <- !is.na(truth[b, a]) && truth[b, a] == 1
        if (!(tab || tba)) next
        cab <- isPos(est[a, b])
        cba <- isPos(est[b, a])
        if (!(cab || cba)) next
        dirTot <- dirTot + 1
        if ((tab && cab && !cba) || (tba && cba && !cab)) dirOk <- dirOk + 1
    }
    pr <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
    re <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
    f1 <- if (!is.na(pr) && !is.na(re) && pr + re > 0)
        2 * pr * re / (pr + re) else NA_real_
    list(tp = tp, fp = fp, fn = fn, tn = tn, unresolved = unresolved,
        precision = pr, recall = re, f1 = f1,
        orientationAccuracy = if (dirTot > 0) dirOk / dirTot else NA_real_,
        coverage = if (nrow(pairs)) resolved / nrow(pairs) else NA_real_)
}
