## Internal engine of causalASCEND. Nothing in this file is exported.
##
## Ancestrality matrix convention (rows and columns are foreground variables):
##   M[a, b] = 1    a is an ancestor of b
##   M[a, b] = 0.5  a is not a descendant of b
##   M[a, b] = 0    neither is an ancestor of the other
##   M[a, b] = NA   unresolved

## Topological order of a 0/1 adjacency matrix (A[i, j] = 1 means i -> j),
## sources first; NULL if the graph has a cycle.
.topoOrder <- function(A) {
    A <- (A != 0) * 1
    A[is.na(A)] <- 0
    k <- nrow(A)
    if (k == 0) return(integer(0))
    indeg <- colSums(A)
    ord <- integer(0)
    zeros <- which(indeg == 0)
    while (length(zeros)) {
        v <- zeros[1]
        zeros <- zeros[-1]
        ord <- c(ord, v)
        nb <- which(A[v, ] != 0)
        for (u in nb) {
            indeg[u] <- indeg[u] - 1
            if (indeg[u] == 0) zeros <- c(zeros, u)
        }
        A[v, nb] <- 0
    }
    if (length(ord) != k) NULL else ord
}

## Write M[i, j] = vij and M[j, i] = vji only if the precedence graph
## (entries 1 or 0.5) stays acyclic.
.tryEdge <- function(M, i, j, vij, vji) {
    Mn <- M
    Mn[i, j] <- vij
    Mn[j, i] <- vji
    b <- (Mn == 1 | Mn == 0.5) * 1
    b[is.na(b)] <- 0
    if (is.null(.topoOrder(b))) list(M = M, ok = FALSE)
    else list(M = Mn, ok = TRUE)
}

## Fisher-z test of the partial correlation of x and y given S, read from
## the correlation matrix C of n samples. Returns a p-value, or NA when the
## test is undefined (too few residual degrees of freedom).
.ciPval <- function(x, y, S, C, n) {
    S <- S[S != x & S != y]
    df <- n - length(S) - 3
    if (df < 1) return(NA_real_)
    m <- C[c(x, y, S), c(x, y, S), drop = FALSE]
    P <- tryCatch(solve(m), error = function(e) tryCatch(
        solve(m + diag(1e-8, nrow(m))), error = function(e2) NULL))
    if (is.null(P)) return(NA_real_)
    den <- sqrt(P[1, 1] * P[2, 2])
    if (!is.finite(den) || den <= 0) return(NA_real_)
    rho <- max(min(-P[1, 2] / den, 1 - 1e-10), -1 + 1e-10)
    2 * stats::pnorm(-abs(sqrt(df) * atanh(rho)))
}

## IAMB Markov-blanket search for target within pool. Because the pool
## holds only non-descendants of target, the blanket is the set of
## nearest ancestors of target within the pool.
.iamb <- function(target, pool, C, n, alpha, maxCond, ci) {
    pool <- setdiff(pool, target)
    mb <- character(0)
    repeat {
        rest <- setdiff(pool, mb)
        if (!length(rest) || length(mb) >= maxCond) break
        pv <- vapply(rest, function(w) ci(target, w, mb, C, n), numeric(1))
        pv[is.na(pv)] <- 1
        k <- which.min(pv)
        if (pv[k] <= alpha) mb <- c(mb, rest[k]) else break
    }
    changed <- TRUE
    while (changed && length(mb)) {
        changed <- FALSE
        for (w in mb) {
            p <- ci(target, w, setdiff(mb, w), C, n)
            if (!is.na(p) && p > alpha) {
                mb <- setdiff(mb, w)
                changed <- TRUE
                break
            }
        }
    }
    mb
}

## Tally R1/R2 (de)activation witnesses for the pair (xi, xj) given the
## guarded set S and return the votes after resolving conflicts.
.witnessVotes <- function(xi, xj, S, C, n, alpha, ci) {
    votes <- c(r1 = 0L, r1r = 0L, r2 = 0L, r2r = 0L)
    evid <- c(r1 = 0, r1r = 0, r2 = 0, r2r = 0)
    ok <- function(p) !is.na(p)
    for (W in S) {
        Sw <- setdiff(S, W)
        pj0 <- ci(W, xj, Sw, C, n, "witness")
        pjI <- ci(W, xj, c(Sw, xi), C, n, "witness")
        pi0 <- ci(W, xi, Sw, C, n, "witness")
        piJ <- ci(W, xi, c(Sw, xj), C, n, "witness")
        if (ok(pj0) && ok(pjI) && pj0 <= alpha && pjI > alpha) {
            votes["r1"] <- votes["r1"] + 1L         # R1: xi precedes xj
            evid["r1"] <- evid["r1"] - log(pj0 + 1e-300)
        }
        if (ok(pi0) && ok(piJ) && pi0 <= alpha && piJ > alpha) {
            votes["r1r"] <- votes["r1r"] + 1L       # R1: xj precedes xi
            evid["r1r"] <- evid["r1r"] - log(pi0 + 1e-300)
        }
        if (ok(pi0) && ok(piJ) && pi0 > alpha && piJ <= alpha) {
            votes["r2"] <- votes["r2"] + 1L         # R2: xi not after xj
            evid["r2"] <- evid["r2"] - log(piJ + 1e-300)
        }
        if (ok(pj0) && ok(pjI) && pj0 > alpha && pjI <= alpha) {
            votes["r2r"] <- votes["r2r"] + 1L       # R2: xj not after xi
            evid["r2r"] <- evid["r2r"] - log(pjI + 1e-300)
        }
    }
    ## opposite directions: keep the one with stronger aggregate evidence
    for (pr in list(c("r1", "r1r"), c("r2", "r2r"))) {
        a <- pr[1]
        b <- pr[2]
        if (votes[a] > 0 && votes[b] > 0) {
            if (evid[a] > evid[b]) votes[b] <- 0L
            else if (evid[b] > evid[a]) votes[a] <- 0L
            else votes[c(a, b)] <- 0L
        }
    }
    if (votes["r1"] > 0) votes["r2"] <- 0L          # strict supersedes weak
    if (votes["r1r"] > 0) votes["r2r"] <- 0L
    votes
}

## The ASCEND algorithm on a correlation matrix C (n samples) whose
## dimnames contain the background labels zlabs and foreground labels
## xlabs. Returns a list with the ancestrality matrix and its provenance.
.ascendEngine <- function(C, n, zlabs, xlabs, maxIter, alpha, alphaMB,
    fdr, minVotes, prescreen, verbose) {
    dX <- length(xlabs)
    vars <- c(zlabs, xlabs)
    isZ <- stats::setNames(vars %in% zlabs, vars)
    maxCond <- max(1, floor((n - 4) / 3))

    ## instrumented CI test: every test of the run goes through ci()
    work <- new.env(parent = emptyenv())
    work$n <- c(pair = 0, witness = 0, mb = 0)
    work$time <- 0
    work$sz <- c(sum = 0, max = 0)
    work$psz <- c(sum = 0, max = 0)
    ci <- function(x, y, S, C, n, kind = "mb") {
        t0 <- Sys.time()
        p <- .ciPval(x, y, S, C, n)
        work$time <- work$time + as.numeric(Sys.time() - t0, units = "secs")
        work$n[kind] <- work$n[kind] + 1
        k <- length(setdiff(S, c(x, y)))
        work$sz <- c(sum = work$sz[["sum"]] + k, max = max(work$sz[["max"]], k))
        if (kind == "pair")
            work$psz <- c(sum = work$psz[["sum"]] + k,
                max = max(work$psz[["max"]], k))
        p
    }

    ## background prescreen: Z marginally associated with each X
    prescreenZ <- function(xi) {
        pv <- vapply(zlabs, function(zc) ci(xi, zc, character(0), C, n),
            numeric(1))
        pv[is.na(pv)] <- 1
        keep <- which(stats::p.adjust(pv, "BH") < prescreen)
        if (length(keep) < min(3, length(zlabs)))
            keep <- order(pv)[seq_len(min(3, length(pv)))]
        zlabs[keep]
    }
    zPool <- stats::setNames(lapply(xlabs, prescreenZ), xlabs)
    nearestAnc <- function(xi, pool)
        .iamb(xi, intersect(pool, vars), C, n, alphaMB, maxCond, ci)
    pa <- stats::setNames(lapply(xlabs, function(xi)
        nearestAnc(xi, zPool[[xi]])), xlabs)

    dn <- list(xlabs, xlabs)
    M <- matrix(NA_real_, dX, dX, dimnames = dn)
    Pp <- matrix(NA_real_, dX, dX, dimnames = dn)
    lastPool <- zPool
    ## provenance of each pair: rule, sweep and number of agreeing witnesses
    prov <- new.env(parent = emptyenv())
    prov$rule <- matrix(NA_character_, dX, dX, dimnames = dn)
    prov$sweep <- matrix(NA_integer_, dX, dX, dimnames = dn)
    prov$votes <- matrix(NA_integer_, dX, dX, dimnames = dn)
    mark <- function(i, j, rule, sweep, v = NA_integer_) {
        prov$rule[i, j] <- prov$rule[j, i] <- rule
        prov$sweep[i, j] <- prov$sweep[j, i] <- as.integer(sweep)
        prov$votes[i, j] <- prov$votes[j, i] <- as.integer(v)
    }

    ## guarded conditioning set: nearest ancestors of either endpoint,
    ## keeping a foreground member only if it is known to precede both
    buildS <- function(i, j, M) {
        cand <- setdiff(union(pa[[xlabs[i]]], pa[[xlabs[j]]]), xlabs[c(i, j)])
        keep <- vapply(cand, function(w) {
            if (isZ[w]) return(TRUE)
            wi <- match(w, xlabs)
            (!is.na(M[wi, i]) && M[wi, i] %in% c(0.5, 1)) &&
                (!is.na(M[wi, j]) && M[wi, j] %in% c(0.5, 1))
        }, logical(1))
        cand[keep]
    }
    symmetry <- function(M, iter, flag) {
        for (i in seq_len(dX)) for (j in seq_len(dX)) {
            if (i != j && !is.na(M[i, j]) && !is.na(M[j, i]) &&
                M[i, j] == 0.5 && M[j, i] == 0.5) {
                M[i, j] <- 0
                M[j, i] <- 0
                mark(i, j, "symmetry", iter)
                flag$changed <- TRUE
            }
        }
        M
    }

    lowerPairs <- which(lower.tri(M), arr.ind = TRUE)
    lowerPairs <- lowerPairs[order(lowerPairs[, "row"], lowerPairs[, "col"]),
        , drop = FALSE]
    converged <- FALSE
    iter <- 0
    while (!converged && iter < maxIter) {
        iter <- iter + 1
        converged <- TRUE
        if (verbose)
            message(sprintf("  [ascend] sweep %d | %d pairs unresolved",
                iter, sum(is.na(M[lower.tri(M)]))))

        ## R3 over all unresolved pairs, one BH correction per sweep
        info <- list()
        for (r in seq_len(nrow(lowerPairs))) {
            i <- lowerPairs[r, "row"]
            j <- lowerPairs[r, "col"]
            if (!is.na(M[i, j])) next
            S <- buildS(i, j, M)
            pv <- ci(xlabs[i], xlabs[j], S, C, n, kind = "pair")
            if (is.na(pv)) {
                converged <- FALSE
                next
            }
            Pp[i, j] <- Pp[j, i] <- pv
            info[[length(info) + 1]] <- list(i = i, j = j, S = S, pv = pv)
        }
        raw <- vapply(info, function(e) e$pv, numeric(1))
        adj <- if (fdr) stats::p.adjust(raw, "BH") else raw
        for (idx in seq_along(info)) {
            e <- info[[idx]]
            i <- e$i
            j <- e$j
            S <- e$S
            if (adj[idx] > alpha) {                   # R3: unrelated
                M[i, j] <- 0
                M[j, i] <- 0
                mark(i, j, "R3", iter)
                converged <- FALSE
                next
            }
            if (!length(S)) next
            votes <- .witnessVotes(xlabs[i], xlabs[j], S, C, n, alpha, ci)
            if (max(votes) < minVotes) next
            res <- NULL
            rule <- NA_character_
            nv <- NA_integer_
            if (votes["r1"] >= minVotes && votes["r1"] >= votes["r1r"]) {
                res <- .tryEdge(M, i, j, 1, 0)
                rule <- "R1"
                nv <- votes[["r1"]]
            }
            if ((is.null(res) || !res$ok) && votes["r1r"] >= minVotes) {
                res <- .tryEdge(M, i, j, 0, 1)
                rule <- "R1"
                nv <- votes[["r1r"]]
            }
            if ((is.null(res) || !res$ok) && votes["r2"] >= minVotes &&
                votes["r2"] >= votes["r2r"]) {
                res <- .tryEdge(M, i, j, 0.5, 0)
                rule <- "R2"
                nv <- votes[["r2"]]
            }
            if ((is.null(res) || !res$ok) && votes["r2r"] >= minVotes) {
                res <- .tryEdge(M, i, j, 0, 0.5)
                rule <- "R2"
                nv <- votes[["r2r"]]
            }
            if (!is.null(res) && res$ok) {
                M <- res$M
                mark(i, j, rule, iter, nv)
                converged <- FALSE
            }
        }

        ## transitive closure: a < k < b implies a < b
        repeat {
            done <- TRUE
            for (k in seq_len(dX)) {
                ancs <- which(M[, k] == 1)
                descs <- which(M[k, ] == 1)
                for (a in ancs) for (b in descs) {
                    if (a != b && (is.na(M[a, b]) || M[a, b] != 1)) {
                        r <- .tryEdge(M, a, b, 1, 0)
                        if (r$ok) {
                            M <- r$M
                            mark(a, b, "transitivity", iter)
                            done <- FALSE
                            converged <- FALSE
                        }
                    }
                }
            }
            if (done) break
        }

        ## a not after b and b not after a implies unrelated
        flag <- new.env(parent = emptyenv())
        flag$changed <- FALSE
        M <- symmetry(M, iter, flag)
        if (flag$changed) converged <- FALSE

        ## refresh nearest ancestors over the enlarged non-descendant pools;
        ## IAMB is deterministic, so skip variables whose pool is unchanged
        for (i in seq_len(dX)) {
            xi <- xlabs[i]
            old <- pa[[xi]]
            pool <- setdiff(unique(c(zPool[[xi]],
                xlabs[which(M[, i] %in% c(0.5, 1))])), xi)
            if (!is.null(lastPool[[xi]]) && setequal(pool, lastPool[[xi]]))
                next
            lastPool[[xi]] <- pool
            pa[[xi]] <- nearestAnc(xi, pool)
            if (!setequal(pa[[xi]], old)) converged <- FALSE
        }
    }

    ## final R3 pass on pairs still unresolved
    for (r in seq_len(nrow(lowerPairs))) {
        i <- lowerPairs[r, "row"]
        j <- lowerPairs[r, "col"]
        if (!is.na(M[i, j])) next
        pv <- ci(xlabs[i], xlabs[j], buildS(i, j, M), C, n, kind = "pair")
        if (!is.na(pv)) Pp[i, j] <- Pp[j, i] <- pv
        if (!is.na(pv) && pv > alpha) {
            M[i, j] <- 0
            M[j, i] <- 0
            mark(i, j, "R3_final", iter)
        }
    }
    M <- symmetry(M, iter, new.env(parent = emptyenv()))

    ## list variables in a topological order of the claimed precedences
    sup <- (!is.na(M) & M %in% c(0.5, 1)) * 1
    dim(sup) <- dim(M)
    diag(sup) <- 0
    ord <- .topoOrder(sup)
    if (!is.null(ord)) M <- M[ord, ord, drop = FALSE]
    diag(M) <- NA
    if (verbose) message(sprintf("  [ascend] converged after %d sweep(s)",
        iter))
    o <- rownames(M)
    nTests <- sum(work$n)
    list(
        M = M,
        pairP = Pp[o, o, drop = FALSE],
        rule = prov$rule[o, o, drop = FALSE],
        ruleSweep = prov$sweep[o, o, drop = FALSE],
        ruleVotes = prov$votes[o, o, drop = FALSE],
        stats = list(
            nTests = nTests,
            nPairTests = unname(work$n["pair"]),
            nWitnessTests = unname(work$n["witness"]),
            nBlanketTests = unname(work$n["mb"]),
            testTimeSec = work$time,
            meanCondSize = if (nTests) work$sz[["sum"]] / nTests else 0,
            maxCondSize = work$sz[["max"]],
            meanPairCondSize = if (work$n[["pair"]])
                work$psz[["sum"]] / work$n[["pair"]] else 0,
            maxPairCondSize = work$psz[["max"]],
            nSweeps = iter))
}
