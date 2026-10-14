set.seed(123)
sim <- simulateTwoTier(n = 1000, nBackground = 20, nForeground = 8,
    r2 = 0.7, sparsity = 0.7, pCross = 0.1, effect = 1.5)

test_that("ascend() returns a valid AscendResult", {
    res <- ascend(sim$background, sim$foreground)
    expect_s4_class(res, "AscendResult")
    expect_true(validObject(res))
    M <- ancestralMatrix(res)
    expect_setequal(rownames(M), rownames(sim$foreground))
    expect_true(all(is.na(diag(M))))
    expect_true(all(is.na(M) | M %in% c(0, 0.5, 1)))
    st <- ascendStats(res)
    expect_equal(st$nTests,
        st$nPairTests + st$nWitnessTests + st$nBlanketTests)
    expect_output(show(res), "AscendResult")
})

test_that("ascend() recovers most ancestral relations in a sparse graph", {
    set.seed(1)
    easy <- simulateTwoTier(n = 2000, nBackground = 20, nForeground = 15,
        r2 = 0.7, sparsity = 0.9, pCross = 0.2, effect = 0.9)
    ev <- evaluateAncestral(ascend(easy$background, easy$foreground),
        trueAncestral(easy))
    expect_gt(ev$precision, 0.6)
    expect_gt(ev$recall, 0.8)
    expect_gt(ev$orientationAccuracy, 0.8)
    expect_gt(ev$coverage, 0.9)
})

test_that("sample order, data frames and transposition are handled", {
    res <- ascend(sim$background, sim$foreground)
    perm <- sample(ncol(sim$foreground))
    res2 <- ascend(sim$background, sim$foreground[, perm])
    expect_identical(ancestralMatrix(res), ancestralMatrix(res2))
    res3 <- ascend(as.data.frame(sim$background),
        as.data.frame(sim$foreground))
    expect_identical(ancestralMatrix(res), ancestralMatrix(res3))
})

test_that("unnamed features get default names", {
    b <- unname(sim$background)
    f <- unname(sim$foreground)
    res <- ascend(b, f)
    expect_setequal(rownames(ancestralMatrix(res)), paste0("x", 1:8))
})

test_that("bad input is rejected", {
    expect_error(ascend(sim$background, sim$foreground[1, , drop = FALSE]),
        "at least two")
    expect_error(ascend(sim$background, sim$foreground[, 1:5]), "same")
    expect_error(ascend(sim$background[, 1:5], sim$foreground[, 1:5]),
        "10 samples")
    expect_error(ascend(sim$background, sim$foreground, alpha = 2), "alpha")
    expect_error(ascend(sim$background, sim$foreground, fdr = NA), "fdr")
    f <- sim$foreground
    rownames(f)[1] <- "z1"
    expect_error(ascend(sim$background, f), "differ")
    expect_error(ascend(letters, sim$foreground), "numeric matrix")
})

test_that("constant features trigger a warning, not an error", {
    f <- sim$foreground
    f[2, ] <- 1
    expect_warning(res <- ascend(sim$background, f), "constant")
    expect_s4_class(res, "AscendResult")
})

test_that("ancestralEdges() lists each relation once", {
    res <- ascend(sim$background, sim$foreground)
    all <- ancestralEdges(res, include = c("ancestor", "non-descendant",
        "unrelated", "unresolved"))
    expect_equal(nrow(all), choose(8, 2))
    expect_named(all, c("from", "to", "relation", "rule", "sweep", "votes",
        "pValue"))
    anc <- ancestralEdges(res, include = "ancestor")
    M <- ancestralMatrix(res)
    expect_equal(nrow(anc), sum(M == 1, na.rm = TRUE))
})

test_that("SummarizedExperiment input works", {
    skip_if_not_installed("SummarizedExperiment")
    se_b <- SummarizedExperiment::SummarizedExperiment(
        assays = list(geno = sim$background))
    se_f <- SummarizedExperiment::SummarizedExperiment(
        assays = list(expr = sim$foreground))
    res <- ascend(se_b, se_f)
    ref <- ascend(sim$background, sim$foreground)
    expect_identical(ancestralMatrix(res), ancestralMatrix(ref))
})

test_that("MultiAssayExperiment input works", {
    skip_if_not_installed("MultiAssayExperiment")
    mae <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = list(snp = sim$background, rna = sim$foreground))
    res <- ascendMultiAssay(mae, background = "snp", foreground = "rna")
    ref <- ascend(sim$background, sim$foreground)
    expect_identical(ancestralMatrix(res), ancestralMatrix(ref))
})
