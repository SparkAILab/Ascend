## The package must return exactly what the implementation used for the
## paper returned. The fixture holds that implementation's output on four
## simulated data sets (generated with the paper's code, R/ascend.R of
## https://github.com/SparkAILab/Ascend).
test_that("ascend() reproduces the reference implementation", {
    ref <- readRDS(test_path("fixtures", "reference_results.rds"))
    for (r in ref) {
        a <- r$args
        set.seed(a$seed)
        sim <- simulateTwoTier(a$n, a$dz, a$dx, r2 = a$r2, sparsity = a$sp,
            pCross = a$pc, effect = a$xe)
        res <- ascend(sim$background, sim$foreground)
        expect_identical(ancestralMatrix(res), r$M)
        expect_equal(pairPValues(res), r$pairP)
        expect_identical(res@rule, r$rule)
    }
})
