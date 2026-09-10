test_that("preprocessCLAMP filters and returns correct structure", {
    Y <- matrix(rnorm(100), 10, 10)
    out <- preprocessCLAMP(Y = Y, mean_cutoff = 0, var_cutoff = 0)
    expect_true(is.list(out))
    expect_true(is.matrix(out$Y_filtered))
    expect_true(is.data.frame(out$rowStats))
})

test_that("preprocessCLAMP errors on non-numeric input", {
    expect_error(preprocessCLAMP(
        Y = data.frame(a = letters[1:5]),
        mean_cutoff = 0, var_cutoff = 0
    ))
})

test_that("preprocessCLAMP log-transforms and fills missing values", {
    Y <- matrix(
        c(0, 3, NA, 100, 127, 255),
        nrow = 2, byrow = TRUE,
        dimnames = list(c("gene1", "gene2"), paste0("sample", 1:3))
    )
    expected <- log2(Y + 1)
    expected[is.na(expected)] <- 0

    expect_message(
        out <- preprocessCLAMP(Y, mean_cutoff = 0, var_cutoff = 0),
        "Applying log2 transformation"
    )
    expect_equal(out$Y_filtered, expected)
    expect_equal(out$rowStats$mean, unname(rowMeans(expected)))
    expect_equal(
        out$rowStats$variance,
        unname(apply(expected, 1, stats::var)) *
            (ncol(expected) - 1) / ncol(expected)
    )
})

test_that("preprocessCLAMP leaves log-scale input unchanged", {
    Y <- matrix(seq(0, 5), nrow = 2)

    expect_message(
        out <- preprocessCLAMP(Y, mean_cutoff = 0, var_cutoff = 0),
        "Already on log scale"
    )
    expect_equal(out$Y_filtered, Y)
})

test_that("preprocessCLAMP log-transforms at the threshold", {
    Y <- matrix(c(0, 99, 1, 100), nrow = 2)

    expect_message(
        out <- preprocessCLAMP(Y, mean_cutoff = 0, var_cutoff = 0),
        "Applying log2 transformation"
    )
    expect_equal(out$Y_filtered, log2(Y + 1))
})

test_that("preprocessCLAMP can skip log2 while cleaning and filtering", {
    Y <- matrix(c(0, 3, NA, 100, 127, 255), nrow = 2, byrow = TRUE)
    expected <- Y
    expected[is.na(expected)] <- 0

    out <- preprocessCLAMP(Y, mean_cutoff = 50, log2_transform = FALSE)
    expect_equal(out$kept_rows, 2L)
    expect_equal(out$Y_filtered, expected[2, , drop = FALSE])
    expect_equal(out$rowStats$mean, mean(expected[2, ]))
    expect_equal(out$rowStats$variance, var(expected[2, ]) * 2 / 3)

    out <- preprocessCLAMP(Y, log2_transform = FALSE)
    expect_equal(out$Y_filtered, expected)
})

test_that("preprocessing rejects invalid log2 flags", {
    Y <- matrix(1:6, nrow = 2)
    for (flag in list(NA, NULL, 1, "FALSE", c(TRUE, FALSE))) {
        expect_error(preprocessCLAMP(Y, log2_transform = flag), "log2_transform")
        expect_error(preprocessCLAMPFBM(Y, log2_transform = flag), "log2_transform")
        expect_error(cleanFBM(Y, log2_transform = flag), "log2_transform")
    }
})


test_that("preprocessCLAMP uses population variance for filtering", {
    Y <- rbind(c(1, 2, 3), c(1, 3, 5), c(2, 2, 2))
    out <- preprocessCLAMP(Y, log2_transform = FALSE)
    expect_equal(out$rowStats$variance, c(2 / 3, 8 / 3, 0))

    filtered <- preprocessCLAMP(Y, var_cutoff = 0.8, log2_transform = FALSE)
    expect_equal(filtered$kept_rows, 2L)
})

test_that("preprocessCLAMP gives zero population variance for one sample", {
    out <- preprocessCLAMP(matrix(c(1, 2), ncol = 1), log2_transform = FALSE)
    expect_equal(out$rowStats$variance, c(0, 0))
    expect_equal(out$kept_rows, 1:2)
})
