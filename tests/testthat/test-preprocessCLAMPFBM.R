test_that("preprocessCLAMPFBM filters FBM and returns correct structure", {
    mat <- matrix(rnorm(100), 10, 10)
    fbm <- FBM(nrow(mat), ncol(mat))
    fbm[, ] <- mat
    out <- preprocessCLAMPFBM(fbm = fbm, mean_cutoff = 0, var_cutoff = 0)
    expect_type(out, "list")
    expect_s4_class(out$fbm_filtered, "FBM")
    expect_setequal(names(out$rowStats), c("row_means", "row_variances"))
    expect_true(is.numeric(out$kept_rows))
})

test_that("preprocessCLAMPFBM errors on non-FBM input", {
    expect_error(preprocessCLAMPFBM(fbm = matrix(1:4, 2, 2)))
})

test_that("FBM preprocessing matches matrix preprocessing with either log2 flag", {
    mat <- matrix(c(0, 3, NA, 100, 127, 255), nrow = 2, byrow = TRUE)
    for (flag in c(TRUE, FALSE)) {
        fbm <- FBM(nrow(mat), ncol(mat), init = mat)
        expected <- preprocessCLAMP(mat, log2_transform = flag)
        out <- preprocessCLAMPFBM(fbm, log2_transform = flag)
        expect_equal(out$fbm_filtered[, ], expected$Y_filtered)
        expect_equal(out$rowStats$row_means, expected$rowStats$mean)
        expect_equal(
            out$rowStats$row_variances,
            expected$rowStats$variance
        )
        expect_equal(out$kept_rows, expected$kept_rows)
        expect_equal(fbm[, ], mat)
    }
})

test_that("cleanFBM defaults to automatic log2 and allows skipping it", {
    mat <- matrix(c(0, 3, NA, 100, 127, 255), nrow = 2)
    for (flag in c(TRUE, FALSE)) {
        fbm <- FBM(nrow(mat), ncol(mat), init = mat)
        if (flag) {
            cleanFBM(fbm)
        } else {
            cleanFBM(fbm, log2_transform = FALSE)
        }
        expected <- if (flag) log2(mat + 1) else mat
        expected[is.na(expected)] <- 0
        expect_equal(fbm[, ], expected)
    }
})
