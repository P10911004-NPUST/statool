#' Title
#'
#' @param data A data frame in which the variables specified in the formula will be found.
#' @param formula A formula specifying the model.
#' @param alpha Numeric value range from 0 to 1 (default: 0.05). The error tolerance.
#' @param p_adjust_method Character string (default: "none"). Other options: `stats::p.adjust.methods`.
#' @param use_art Whether using aligned rank transform. Default: TRUE
#'
#' @return A list, contains omnibus test and post-hoc test results.
#' @import ARTool stats utils
#' @export
oneway_test <- function(
        data,
        formula,
        alpha = 0.05,
        p_adjust_method = "holm",
        use_art = TRUE
){
    p_adjust_method <- match.arg(p_adjust_method, stats::p.adjust.methods)

    df0 <- stats::model.frame(formula, data, drop.unused.levels = TRUE)
    y <- colnames(df0)[1]
    x <- colnames(df0)[2]
    colnames(df0) <- c("y", "x")

    is_normal <- is_normality(df0, y ~ x)  # from "./utils.R"
    # is_normal <- with(df0, is_normality(df0, get(y) ~ get(x)))  # from "./utils.R"
    is_balance <- !is_unbalance(df0, y ~ x)  # from "./utils.R"
    is_var_equal <- levene_test(df0, y ~ x)[["is_var_equal"]]  # from "./homoscedasticity.R"

    ##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ## Parametric ====
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
                alpha = alpha,
                p_adjust_method = p_adjust_method
            )
        }

        if (is_var_equal & !is_balance)
        {
            tests <- "Fisher's ANOVA + Tukey-Kramer test"
            post_hoc <- Tukey_Kramer_test(
                data = df0,
                formula = y ~ x,
                alpha = alpha,
                p_adjust_method = p_adjust_method
            )
        }

        if ( ! is_var_equal )
        {
            tests <- "Welch's ANOVA + Games-Howell test"
            post_hoc <- Games_Howell_test(
                data = df0,
                formula = y ~ x,
                alpha = alpha,
                p_adjust_method = p_adjust_method
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
            post_hoc <- list()

            if ( ! is.factor(df0$x) ) df0$x <- as.factor(df0$x)

            art_mod <- ARTool::art(formula = y ~ x, data = df0)
            pre_hoc <- stats::anova(art_mod)
            pre_hoc_pass <- pre_hoc$`Pr(>F)` < alpha

            art_c <- ARTool::art.con(
                m = art_mod,
                formula = ~ x,
                adjust = p_adjust_method
            )
            art_c <- as.data.frame(art_c)

            #### Tidy-up contrasts
            desc_df <- with(
                data = df0,
                expr = vapply(
                    X = c("length", "mean", "sd", "min", "max", "median"),
                    FUN = function(fns) tapply(y, x, fns),
                    FUN.VALUE = numeric(length(unique(df0[["x"]])))
                )
            )
            desc_df <- as.data.frame(desc_df)
            desc_df <- desc_df[order(desc_df[["median"]], decreasing = TRUE), ]
            group_comb <- utils::combn(rownames(desc_df), m = 2, simplify = FALSE)

            cont_vct <- vector("character", length = length(group_comb))
            pval_vct <- vector("numeric", length = length(group_comb))
            for (i in seq_along(group_comb)) {
                gname1 <- group_comb[[i]][1]
                gname2 <- group_comb[[i]][2]
                bool1 <- grepl(gname1, art_c$contrast)
                bool2 <- grepl(gname2, art_c$contrast)
                pval <- art_c[["p.value"]][bool1 & bool2]

                cont_vct[i] <- paste(gname1, gname2, sep = " |vs| ")
                pval_vct[i] <- pval
            }

            cld_res <- compact_letter_display(
                groups = row.names(desc_df),
                means = desc_df[["median"]],
                comparisons = cont_vct,
                pval = pval_vct,
                alpha = alpha
            )

            result <- data.frame(
                row.names = row.names(desc_df),
                GROUPS = row.names(desc_df),
                N = desc_df[["length"]],
                AVG = desc_df[["mean"]],
                SD = desc_df[["sd"]],
                MED = desc_df[["median"]],
                MIN = desc_df[["min"]],
                MAX = desc_df[["max"]],
                CLD = unname(cld_res)
            )

            post_hoc$result <- result
            post_hoc$comparisons <- art_c
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
            pre_hoc_pass <- pre_hoc[["p.value"]] < alpha

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




