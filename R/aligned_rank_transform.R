#' One-way aligned rank transform ANOVA
#'
#' @param data A data frame.
#' @param formula A formula specifying the model, for example, yield ~ fertilizer_level.
#'
#' @importFrom car Anova
#'
#' @return ANOVA (type 3) result
#' @export
#'
#' @examples
#' set.seed(1)
#' df0 <- data.frame(
#'     group = as.factor(rep(c("A", "B", "C"), each = 10)),
#'     norm_data = c(rnorm(10, -1, 2), rnorm(10, 3, 2), rnorm(10, 0, 1.5)),
#'     skew_data = c(sqrt(rnorm(10, -7.2, 2) ^ 2), runif(10), runif(10))
#' )
#' art_1(df0, skew_data ~ group)
#' #> Anova Table (Type 3)
#' #> Response: group
#' #> Sum Sq Df F value    Pr(>F)
#' #> (Intercept) 6502.5  1 263.022 1.917e-15 ***
#' #> x           1580.0  2  31.955 7.623e-08 ***
#' #> Residuals    667.5 27
#' #> ---
#' #> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
art_1 <- function(data, formula)
{
    df0 <- stats::model.frame(formula, data)
    colnames(df0) <- c("y", "x")
    df0 <- df0[!is.na(df0$y), ]

    grand_mean <- mean(df0$y)

    aov_mod <- stats::aov(y ~ x, df0)
    df0$residual <- stats::residuals(aov_mod)

    group_means <- with(df0, tapply(y, x, mean))
    estimated_effect <- group_means - grand_mean

    df0$estimated_effect <- vapply(
        X = df0$x,
        FUN = function(x) estimated_effect[match(x, names(estimated_effect))],
        FUN.VALUE = numeric(1)
    )

    df0$y_prime <- with(df0, residual + estimated_effect)  # aligned response
    df0$`Y''` <- rank(df0$y_prime, ties.method = "average")  # ranked response

    aov_mod <- stats::aov(`Y''` ~ x, df0)
    pre_hoc <- car::Anova(aov_mod, type = 3)
    attr(pre_hoc, "heading") <- sprintf("Anova Table (Type 3)\nResponse: %s", as.character(formula)[3])
    attr(pre_hoc, "data") <- df0
    attr(pre_hoc, "estimated_effect") <- df0$estimated_effect

    return(pre_hoc)

    ##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ## Test
    ##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # set.seed(1)
    #
    # df0 <- data.frame(
    #     group = as.factor(rep(c("A", "B", "C"), each = 10)),
    #     norm_data = c(rnorm(10, -1, 2), rnorm(10, 3, 2), rnorm(10, 0, 1.5)),
    #     skew_data = c(sqrt(rnorm(10, -7.2, 2) ^ 2), runif(10), runif(10))
    # )
    #
    # aov_mod <- aov(skew_data ~ group, df0)
    # car::Anova(aov_mod, type = 3)
    # art_mod <- ARTool::art(skew_data ~ group, df0)
    # stats::anova(art_mod)
    #
    # mod <- art_1(df0, skew_data ~ group)
    # attr(mod, "data")
    #
    # em <- emmeans::emmeans(aov(skew_data ~ group, df0), "group")
    # em
}




