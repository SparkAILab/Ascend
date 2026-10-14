#' Simulate two-tier data
#'
#' Draws data from a linear-Gaussian structural equation model with a
#' background layer that precedes a foreground layer: background
#' variables form a random DAG, each foreground variable receives edges
#' from random background variables, and foreground variables form a
#' random DAG among themselves. No foreground variable affects a
#' background variable.
#'
#' @param n Number of samples.
#' @param nBackground,nForeground Numbers of background and foreground
#'   variables.
#' @param r2 Fraction of each foreground variable's variance explained by
#'   its parents.
#' @param sparsity Sparsity of the within-layer DAGs (larger is sparser;
#'   an edge appears with probability `1 - sparsity`, kept in
#'   `[0.01, 0.4]`).
#' @param pCross Probability of each background-to-foreground edge.
#' @param effect Size of the foreground-to-foreground effects, relative to
#'   the parent's standard deviation.
#' @param linear Fraction of parent effects that are linear; the rest are
#'   squared, square-root, softplus or ReLU transforms of the parent.
#' @param scaleBackground Logical: standardise each background variable
#'   to unit variance as it is drawn. Without it, background variances
#'   grow along background paths and swamp the foreground signal when the
#'   background is large.
#'
#' @details Call [set.seed()] first for reproducible data.
#'
#' @return A list with `background` and `foreground` (numeric matrices,
#'   features in rows, samples in columns), `adjacency` (the true
#'   foreground graph, `adjacency[parent, child] = 1`) and `params`.
#' @seealso [trueAncestral()], [ascend()]
#' @examples
#' set.seed(1)
#' sim <- simulateTwoTier(n = 200, nBackground = 5, nForeground = 4)
#' dim(sim$background)
#' sim$adjacency
#' @export
simulateTwoTier <- function(n, nBackground, nForeground, r2 = 0.5,
    sparsity = 0.9, pCross = 0.05, effect = 0.8,
    linear = 1, scaleBackground = TRUE) {
    .checkScalar(n, "n", 2, Inf)
    .checkScalar(nBackground, "nBackground", 1, Inf)
    .checkScalar(nForeground, "nForeground", 1, Inf)
    .checkScalar(r2, "r2", 0, 1)
    .checkScalar(sparsity, "sparsity", 0, 1)
    .checkScalar(pCross, "pCross", 0, 1)
    .checkScalar(linear, "linear", 0, 1)
    dZ <- nBackground
    dX <- nForeground

    prep <- function(mat, pr) {
        mat <- as.matrix(mat)
        p <- ncol(mat)
        if (p == 0 || pr >= 1) return(mat)
        nNl <- round((1 - pr) * p)
        if (nNl <= 0) return(mat)
        cols <- sample.int(p, nNl)
        types <- sample(c("sq", "sqrt", "softplus", "relu"), nNl,
            replace = TRUE)
        for (ii in seq_along(cols)) {
            v <- mat[, cols[ii]]
            mat[, cols[ii]] <- switch(types[ii], sq = v^2,
                sqrt = sqrt(abs(v)), softplus = log1p(exp(v)),
                relu = pmax(v, 0))
        }
        mat
    }
    noiseSd <- function(signal, r2) {
        v <- stats::var(signal)
        if (is.na(v) || v == 0) return(1)
        if (r2 <= 0) return(sqrt(v))
        sqrt(max(0, (v - r2 * v) / r2))
    }
    randDag <- function(k) {
        pEdge <- min(max(1 - sparsity, 0.01), 0.4)
        A <- matrix(0, k, k)
        if (k > 1) for (i in seq_len(k - 1)) for (j in (i + 1):k)
            if (stats::runif(1) < pEdge) A[i, j] <- 1
        perm <- sample.int(k)
        A[perm, perm]
    }

    Azz <- randDag(dZ)
    Axx <- randDag(dX)
    Azx <- matrix(stats::rbinom(dZ * dX, 1, pCross), dZ, dX)

    Z <- matrix(0, n, dZ)
    for (idx in .topoOrder(Azz)) {
        par <- which(Azz[, idx] == 1)
        Z[, idx] <- if (!length(par)) stats::rnorm(n) else
            as.numeric(prep(Z[, par, drop = FALSE], linear) %*%
                stats::rnorm(length(par))) + stats::rnorm(n)
        if (scaleBackground)
            Z[, idx] <- (Z[, idx] - mean(Z[, idx])) / stats::sd(Z[, idx])
    }

    X <- matrix(0, n, dX)
    adj <- matrix(0, dX, dX)
    ordX <- .topoOrder(Axx)
    for (s in seq_along(ordX)) {
        j <- ordX[s]
        zp <- which(Azx[, j] == 1)
        sigZ <- if (length(zp))
            as.numeric(prep(Z[, zp, drop = FALSE], linear) %*%
                sample(c(1, -1), length(zp), TRUE)) else rep(0, n)
        xp <- intersect(which(Axx[, j] == 1),
            if (s > 1) ordX[seq_len(s - 1)] else integer(0))
        sigX <- rep(0, n)
        if (length(xp)) {
            Pax <- prep(X[, xp, drop = FALSE], linear)
            adj[xp, j] <- 1
            csd <- apply(Pax, 2, stats::sd)
            csd[csd == 0] <- 1
            sigX <- as.numeric(Pax %*%
                ((effect / csd) * sample(c(1, -1), length(xp), TRUE)))
        }
        sig <- sigZ + sigX
        if (stats::var(sig) < 1e-6) sig <- stats::rnorm(n)
        X[, j] <- sig + stats::rnorm(n, sd = noiseSd(sig, r2))
    }

    zn <- paste0("z", seq_len(dZ))
    xn <- paste0("x", seq_len(dX))
    sn <- paste0("sample", seq_len(n))
    background <- t(Z)
    foreground <- t(X)
    dimnames(background) <- list(zn, sn)
    dimnames(foreground) <- list(xn, sn)
    dimnames(adj) <- list(xn, xn)
    list(background = background, foreground = foreground, adjacency = adj,
        params = list(n = n, nBackground = dZ, nForeground = dX, r2 = r2,
            sparsity = sparsity, pCross = pCross, effect = effect,
            linear = linear, scaleBackground = scaleBackground))
}
