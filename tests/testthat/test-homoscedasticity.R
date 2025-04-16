test_that(
    desc = "Levene test for testing data homoscedasticity",
    code = {
        utils::data("iris", package = "datasets")
        out <- levene_test(iris, Petal.Length ~ Species)

        testthat::expect_equal(out$is_var_equal, FALSE)
    }
)
