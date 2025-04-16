test_that(
    desc = "multiplication works",
    code = {
        set.seed(1)
        df0 <- data.frame(
            group = as.factor(rep(c("A", "B", "C"), each = 10)),
            skew_data = c(rlnorm(10, -1, 2), sqrt(rcauchy(10, 3, 2) ^ 2), rlnorm(10, 0, 1.5))
        )
        out <- art_1(df0, skew_data ~ group)

        testthat::expect_equal(unname(round(out$pre_hoc$statistic, 4)), 6.851)
        testthat::expect_equal(round(out$pre_hoc$p.value, 4), 0.0039)
        testthat::expect_equal(out$cld$CLD, c("a", "ab", "b"))
    }
)
