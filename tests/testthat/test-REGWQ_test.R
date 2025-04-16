testthat::test_that(
    desc = "REGWQ test for normal distributed, balance designed, and homoscedastic data.",
    code = {
        set.seed(1)
        group <- rep(c("A", "B", "C", "D", "E"), each = 10)
        val <- c(rnorm(10, 1, 1.5), rnorm(10, 2, 1), rnorm(10, -4, 2), rnorm(10), rnorm(10))
        df0 <- data.frame(group = group, val = val)
        out <- REGWQ_test(df0, val ~ group)

        # Check group names
        testthat::expect_identical(out$cld$GROUP, c("B", "A", "E", "D", "C"))

        # Check CLD
        testthat::expect_identical(out$cld$CLD, c("a", "ab", "b", "b", "c"))

        # Check padj
        testthat::expect_equal(
            object = round(out$post_hoc$padj, 4),
            expected = c(0.0561, 0.0008, 0.0531, 0.0014, 0.1211, 0.9801, 0, 0, 0, 0)
        )
    }
)



