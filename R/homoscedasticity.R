# Test if the data variances are equal (homoscedasticity)

# For grouped data ====
## F-test
## Barlett's test
## Levene's test
## Brown-Forsythe test (similar with Levene's test, but using median)
## Fligner-Killeen test
## Cochran's C test
## Conover's test
## Hartley's test
## O'Brien's test


#' Homoscedasticity
#'
#' @description
#' Test if the variances among groups are equal.
#'
#' @param data data frame
#' @param formula formula
#' @param alpha (default: 0.05)
#' @param method Available options are: "median" (default), "mean", and "trim_mean".
#'
#' @returns Logical value
#' @export
#'
#' @examples
#' utils::data("iris", package = "datasets")
#' is_var_equal(iris, Petal.Length ~ Species)
#' #> [1] FALSE
is_var_equal <- function(data, formula, alpha = 0.05, method = "median")
{
    method <- match.arg(method, c("median", "mean", "trim_mean", "F"))

    if (method == "F")
    {
        ret <- F_test(data, formula)
        return(ret[["pval"]] > alpha)
    }

    levene_test(data, formula, alpha, method)[["is_var_equal"]]
}

#' @title Levene's test
#' @description Computes Levene's test for homogeneity of variance across groups.
#'
#' @param data data frame
#' @param formula formula
#' @param alpha (default: 0.05)
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
levene_test <- function(data, formula, alpha = 0.05, method = "median")
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
    is_var_equal <- pval > alpha

    ret <- list(
        "is_var_equal" = is_var_equal,
        "statistics" = aov_mod,
        "pval" = pval
    )

    return(ret)
}


F_test <- function(
        data,
        formula = NULL,
        x = NULL,
        y = NULL,
        alternative = c("two.sided", "less", "greater"),
        alpha = 0.05
) {
    alternative <- match.arg(alternative, c("two.sided", "less", "greater"))

    if (is.null(formula))
    {
        if (inherits(x, "character"))
            x <- data[[x]]

        if (inherits(y, "character"))
            y <- data[[y]]

        fstats <- stats::var.test(
            x = x,
            y = y,
            alternative = alternative,
            conf.level = 1 - alpha
        )

    } else {
        # df0 <- stats::model.frame(formula, data, drop.unused.levels = TRUE)
        # x_name <- colnames(df0)[2]
        # y_name <- colnames(df0)[1]
        # colnames(df0) <- c("y", "x")
        # df0 <- df0[stats::complete.cases(df0[["y"]]), ]

        fstats <- stats::var.test(
            formula,
            data,
            alternative = alternative,
            conf.level = 1 - alpha
        )
    }

    fstats <- data.frame(
        comparison = fstats[["data.name"]],
        `F` = unname(fstats[["statistic"]]),
        df_num = fstats[["parameter"]][["num df"]],
        df_denom = fstats[["parameter"]][["denom df"]],
        pval = fstats[["p.value"]],
        conf_int_lower = fstats[["conf.int"]][1],
        conf_int_upper = fstats[["conf.int"]][2],
        var_ratio = unname(fstats[["estimate"]])
    )

    return(fstats)
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
barlett_test <- function(data, formula)
{
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
fligner_test <- function(data, formula)
{
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


