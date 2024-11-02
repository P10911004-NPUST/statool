# Test if the data variances are equal (homoscedasticity)

## There are several solutions to test for the equality of variance across groups, including:
## 1. F-test
## 2. Barlett's test
## 3. Levene's test
## 4. Fligner-Killeen test


F_test <- function(data, formula){
    cat("Not yet")
}


barlett_test <- function(data, formula){
    # If the data input is a dataframe
    if (is.data.frame(data)){
        y_name <- as.character(formula)[2]
        x_name <- as.character(formula)[3]
        # res <- bartlett.test(get(y_name) ~ get(x_name), data)

        df0 <- data.frame(
            y = data[[y_name]],
            x = data[[x_name]]
        )

        df0 <- df0[!is.na(df0$x), ]

        res <- stats::bartlett.test(formula = y ~ x, data = df0)
    } else{
        cat("Data input should be a dataframe")
    }

    return(res$p.value)
}




#' Levene's test
#'
#' Computes Levene's test for homogeneity of variance across groups (homoscedasticity).
#'
#' @param data data frame
#' @param formula formula
#' @param method options are: "trim_mean" (default), "mean", "median"
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
#' levene_test(iris, Petal.Length ~ Species)
#' #> $is_var_equal
#' #> [1] FALSE
#' #> $statistic
#' #> Analysis of Variance Table
#' #>
#' #> Response: diff
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
levene_test <- function(data, formula, method = NULL){
    if (is.null(method)) method <- "trim_mean"
    method <- match.arg(method, c("trim_mean", "mean", "median"))

    trim_mean <- function(x) mean(x, trim = 0.1)

    df0 <- stats::model.frame(formula = formula, data = data)
    colnames(df0) <- c("y", "x")

    xbar <- stats::aggregate(x = df0$y, by = list(df0$x), FUN = get(method))
    colnames(xbar) <- c("x", "xbar")

    df0 <- merge(df0, xbar, by.x = "x")
    df0$diff <- abs(df0$y - df0$xbar)

    aov_mod <- stats::aov(formula = diff ~ x, data = df0)
    aov_mod <- stats::anova(aov_mod)

    pval <- aov_mod$`Pr(>F)`[1]
    is_var_equal <- pval > 0.05

    ret <- list(
        "is_var_equal" = is_var_equal,
        statistic = aov_mod,
        pval = pval
    )

    return(ret)
}




fligner_test <- function(data, formula){
    # If the data input is a dataframe
    if (is.data.frame(data)){
        y_name <- as.character(formula)[2]
        x_name <- as.character(formula)[3]
        # res <- bartlett.test(get(y_name) ~ get(x_name), data)

        df0 <- data.frame(
            y = data[[y_name]],
            x = data[[x_name]]
        )

        df0 <- df0[!is.na(df0$x), ]

        res <- stats::fligner.test(formula = y ~ x, data = df0)
    } else{
        cat("Data input should be a dataframe")
    }

    return(res$p.value)
}


