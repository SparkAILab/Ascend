test_that("simulateTwoTier() is reproducible and respects the tiers", {
    set.seed(9)
    a <- simulateTwoTier(200, 5, 6, sparsity = 0.5)
    set.seed(9)
    b <- simulateTwoTier(200, 5, 6, sparsity = 0.5)
    expect_identical(a, b)
    expect_equal(dim(a$background), c(5, 200))
    expect_equal(dim(a$foreground), c(6, 200))
    expect_equal(unname(apply(a$background, 1, stats::sd)), rep(1, 5))
    expect_true(all(a$adjacency %in% c(0, 1)))
})

test_that("trueAncestral() closes paths transitively", {
    A <- matrix(0, 3, 3, dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
    A["a", "b"] <- A["b", "c"] <- 1
    R <- trueAncestral(A)
    expect_equal(R["a", "c"], 1)
    expect_equal(R["c", "a"], 0)
    expect_true(is.na(R["a", "a"]))
})

test_that("evaluateAncestral() scores a perfect estimate as perfect", {
    set.seed(9)
    sim <- simulateTwoTier(200, 5, 6, sparsity = 0.5)
    truth <- trueAncestral(sim)
    ev <- evaluateAncestral(truth, truth)
    expect_equal(ev$fp, 0)
    expect_equal(ev$fn, 0)
    expect_equal(ev$coverage, 1)
    rev <- t(truth)
    ev2 <- evaluateAncestral(rev, truth)
    if (ev$tp > 0) expect_equal(ev2$orientationAccuracy, 0)
})
