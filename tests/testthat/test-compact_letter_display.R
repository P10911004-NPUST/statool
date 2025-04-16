test_that(
    desc = "multiplication works",
    code = {
        utils::data("iris", package = "datasets")
        avg <- tapply(iris$Sepal.Length, iris$Species, "mean")
        aov_mod <- stats::aov(Sepal.Length ~ Species, iris)
        res <- stats::TukeyHSD(aov_mod)
        res <- as.data.frame(res$Species)

        out <- compact_letter_display(
            groups = names(avg),
            means = avg,
            comparisons = row.names(res),
            pval = res$`p adj`,
            comparison_symbol = "-"
        )

        testthat::expect_equal(names(out), c("setosa", "versicolor", "virginica"))
        testthat::expect_equal(unname(out), c("c", "b", "a"))
    }
)
