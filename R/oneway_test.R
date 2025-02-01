#' Title
#'
#' @param data A data frame in which the variables specified in the formula will be found.
#' @param formula A formula specifying the model.
#' @param alpha Numeric value range from 0 to 1 (default: 0.05). The error tolerance.
#' @param p_adjust_method Character string (default: "none"). Other options: `stats::p.adjust.methods`.
#' @param use_art Whether using aligned rank transform. Default: TRUE
#'
#' @return A list, contains omnibus test and post-hoc test results.
#' @import ARTool
#' @export
oneway_test <- function(
        data,
        formula,
        alpha = 0.05,
        p_adjust_method = "holm",
        use_art = FALSE
){
    p_adjust_method <- base::match.arg(p_adjust_method, stats::p.adjust.methods)

    df0 <- stats::model.frame(formula, data, drop.unused.levels = TRUE)
    y_name <- colnames(df0)[1]
    x_name <- colnames(df0)[2]
    colnames(df0) <- c("y", "x")

    desc_df <- with(
        data = df0,
        expr = vapply(
            X = c("length", "mean", "sd", "min", "max", "median"),
            FUN = function(fns) tapply(y, x, fns),
            FUN.VALUE = numeric(length(unique(df0$x)))
        )
    )

    is_normal <- is_normality(df0, y ~ x)
    is_balance <- !is_unbalance(df0, y ~ x)
    is_var_equal <- levene_test(df0, y ~ x)$is_var_equal

    ##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ## Parametric ====
    ## Note: p-values were not adjusted in parametric test
    ##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if (is_normal)
    {
        pre_hoc <- stats::oneway.test(y ~ x, df0, var.equal = is_var_equal)
        pre_hoc_pass <- pre_hoc$p.value < alpha

        if (is_var_equal & is_balance)
        {
            tests <- "Fisher's ANOVA + Tukey-HSD test"
            post_hoc <- Tukey_HSD_test(
                data = df0,
                formula = y ~ x,
                alpha = alpha
            )
        }

        if (is_var_equal & !is_balance)
        {
            tests <- "Fisher's ANOVA + Tukey-Kramer test"
            post_hoc <- Tukey_Kramer_test(
                data = df0,
                formula = y ~ x,
                alpha = alpha
            )
        }

        if ( ! is_var_equal )
        {
            tests <- "Welch's ANOVA + Games-Howell test"
            post_hoc <- Games_Howell_test(
                data = df0,
                formula = y ~ x,
                alpha = alpha
            )
        }
    }

    ##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ## Non-parametric ====
    ##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if ( ! is_normal )
    {
        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ### Aligned Rank Transformed (ART) + ART-C ====
        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (use_art)
        {
            tests <- "Aligned Rank Transform (ART) + ART-Contrast"

            art_mod <- ARTool::art(formula = y ~ x, data = df0)
            pre_hoc <- stats::anova(art_mod)
            pre_hoc_pass <- pre_hoc$`Pr(>F)` < alpha

            post_hoc <- ARTool::art.con(art_mod, ~ x, adjust = p_adjust_method)
            post_hoc <- as.data.frame(post_hoc)
        }

        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ### Kruskal-Wallis + Dunn's test ====
        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if ( ! use_art )
        {
            tests <- "Kruskal-Wallis + Dunn's test"

            pre_hoc <- with(
                data = df0,
                expr = stats::kruskal.test(y, x)
            )
            pre_hoc_pass <- pre_hoc$p.value < alpha

            post_hoc <- Dunn_test(
                data = df0,
                formula = y ~ x,
                alpha = alpha,
                p_adjust_method = p_adjust_method
            )
        }
    }

    ret <- list(
        pre_hoc = pre_hoc,
        result = post_hoc$result,
        comparison = post_hoc$comparisons,
        tests = tests
    )

    return(ret)
}




