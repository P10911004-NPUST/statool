# Test if the data variances are equal (homoscedasticity)

# For grouped data ====
## F-test
## Barlett's test
## Levene's test (including Brown-Forsythe test, when using median or trimmed-mean)
## Fligner-Killeen test
## Cochran's C test
## Conover's test
## Hartley's test
## O'Brien's test


F_test <- function(data, formula){
    cat("Not yet")
}


#' Barlett's test
#'
#' This is a wrapper for the `stats::bartlett.test()` function.
#' Computes Barlett's test for homogeneity of variance across groups (homoscedasticity).
#'
#' @param data data frame
#' @param formula formula
#'
#' @return returns a list containing three elements:
#' 1. is_var_equal
#' 2. statistic
#' 3. pval
#' @export
#' @author Joon-Keat Lai
#'
#' @import datasets
#' @importFrom utils data
#'
#' @examples
#' utils::data("iris", package = "datasets")
#' barlett_test(iris, Petal.Length ~ Species)
#' #> $is_var_equal
#' #> [1] FALSE
#' #>
#' #> $statistic
#' #> Bartlett's K-squared
#' #>              55.4225
#' #>
#' #> $pval
#' #> [1] 2.005419e-08
#'
#' @seealso [stats::bartlett.test()]
barlett_test <- function(data, formula){
    if (!is.data.frame(data)) stop("`data` should be a dataframe")

    df0 <- stats::model.frame(formula = formula, data = data)
    colnames(df0) <- c("y", "x")

    df0 <- df0[stats::complete.cases(df0$x), , drop = FALSE]

    results <- stats::bartlett.test(formula = y ~ x, data = df0)

    pval <- results$p.value
    is_var_equal <- pval > 0.05
    statistic <- results$statistic

    ret <- list(
        "is_var_equal" = is_var_equal,
        "statistic" = statistic,
        "pval" = pval
    )

    return(ret)
}



#' Levene's test
#'
#' Computes Levene's test for homogeneity of variance across groups (homoscedasticity).
#'
#' @param data data frame
#' @param formula formula
#' @param method Available options are: "median" (default), "mean", and "trim_mean".
#'
#' @return returns a list containing three elements:
#' 1. is_var_equal
#' 2. statistic
#' 3. pval
#'
#' @export
#'
#' @author Joon-Keat Lai
#'
#' @import datasets
#' @import utils
#'
#' @examples
#' utils::data("iris", package = "datasets")
#' levene_test(iris, Petal.Length ~ Species)
#'
#' #> $is_var_equal
#' #> [1] FALSE
#' #>
#' #> $statistic
#' #> Analysis of Variance Table (trim_mean)
#' #>
#' #> Response: abs(x - trim_mean)
#' #>            Df Sum Sq Mean Sq F value    Pr(>F)
#' #> x           2 2.6526 1.32629  20.045 2.005e-08 ***
#' #> Residuals 147 9.7265 0.06617
#' #> ---
#' #> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#' #>
#' #> $pval
#' #> [1] 2.005419e-08
#'
#' @references
#' Fox J. and Weisberg S. 2019.
#' \emph{An R Companion to Applied Regression},
#' Third Edition, Sage.
#' @seealso [car::leveneTest()]
levene_test <- function(data, formula, method = "median")
{
    method <- match.arg(method, c("median", "mean", "trim_mean"))

    trim_mean <- function(x) mean(x, trim = 0.1)

    df0 <- stats::model.frame(formula = formula, data = data)
    colnames(df0) <- c("y", "x")

    # xbar could be median, mean, or trimmed mean
    xbar <- stats::aggregate(x = df0$y, by = list(df0$x), FUN = get(method))
    colnames(xbar) <- c("x", "xbar")

    df0 <- merge(df0, xbar, by.x = "x")
    df0$diff <- abs(df0$y - df0$xbar)

    aov_mod <- stats::aov(formula = diff ~ x, data = df0)
    aov_mod <- stats::anova(aov_mod)
    attr(aov_mod, "heading") <- c(
        sprintf("Analysis of Variance Table (%s)\n", method),
        sprintf("Response: abs(x - %s)", method)
    )

    pval <- aov_mod$`Pr(>F)`[1]
    is_var_equal <- pval > 0.05

    ret <- list(
        "is_var_equal" = is_var_equal,
        "statistics" = aov_mod,
        "pval" = pval
    )

    return(ret)
}



#' Fligner's test
#'
#' This is a wrapper for the `stats::fligner.test()` function.
#' Computes Fligner's test for homogeneity of variance across groups (homoscedasticity).
#'
#' @param data data frame
#' @param formula formula
#'
#' @return returns a list containing three elements:
#' 1. is_var_equal
#' 2. statistic
#' 3. pval
#' @export
#' @author Joon-Keat Lai
#'
#' @import datasets
#' @importFrom utils data
#'
#' @examples
#' utils::data("iris", package = "datasets")
#' fligner_test(iris, Petal.Length ~ Species)
#' #> $is_var_equal
#' #> [1] FALSE
#' #>
#' #> $statistic
#' #> Fligner-Killeen:med chi-squared
#' #>                        34.70802
#' #>
#' #> $pval
#' #> [1] 2.90569e-08
#'
#' @seealso [stats::fligner.test()]
fligner_test <- function(data, formula){
    if (!is.data.frame(data)) stop("`data` should be a dataframe")

    df0 <- stats::model.frame(formula = formula, data = data)
    colnames(df0) <- c("y", "x")

    df0 <- df0[stats::complete.cases(df0$x), , drop = FALSE]

    results <- stats::fligner.test(formula = y ~ x, data = df0)

    pval <- results$p.value
    is_var_equal <- pval > 0.05
    statistic <- results$statistic

    ret <- list(
        "is_var_equal" = is_var_equal,
        "statistic" = statistic,
        "pval" = pval
    )

    return(ret)
}



# For regression ====
## Goldfeld-Quandt test
## Park test
## Glejser test
## Harrison-McCabe test
## Breusch-Pagan test
## White test
## Cook-Weisberg test


