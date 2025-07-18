#' Games-Howell test
#'
#' Represent significance statements resulting from all-pairwise comparisons (when group variances are unequal).
#'
#' @param data A data frame.
#' @param formula A formula specifying the model, for example, yield ~ fertilizer_level.
#' @param alpha Numeric value range from 0 to 1 (default: 0.05). The error tolerance.
#' @param p_adjust_method A character string (default: "none"). Other options: `stats::p.adjust.methods`.
#'
#' @return A list with four elements:
#' 1. tests: A message showing the statistical methods applied on the dataset.
#' 2. pre_hoc: The result of pmnibus test.
#' 3. post_hoc: includes statistics parameters for each pairwise comparisons.
#' 4. cld: A dataframe reporting the descriptive stats and compact letter display.
#'
#' @export
#' @author Joon-Keat Lai
#'
#' @references
#' Howell, D.C. 2013.
#' Statistical methods for psychology (8th ed.). pg. 393.
#' Wadsworth Cengage Learning, Belmont, CA.
#'
#' @examples
#' set.seed(1)
#' df0 <- data.frame(
#'     group = rep(c("A", "B", "C"), each = 10),
#'     norm_data = c(rnorm(10, -1, 3), rnorm(10, 3, 1), rnorm(10, 0, 1)),
#'     skew_data = c(sqrt(rnorm(10, -7.2, 2) ^ 2), runif(10), runif(10))
#' )
#'
#' out <- Games_Howell_test(df0, norm_data ~ group)
#' out$cld
#' #>   GROUP  N        AVG        SD          MED        MIN       MAX CLD
#' #> 1     B 10  3.2488450 1.0695148  3.491872279  0.7853001 4.5117812   a
#' #> 2     C 10 -0.1336732 0.9556076  0.009218122 -1.9893517 0.9189774   b
#' #> 3     A 10 -0.6033917 2.3417579 -0.230273356 -3.5068858 3.7858424   b
Games_Howell_test <- function(
        data,
        formula,
        alpha = 0.05,
        p_adjust_method = "none"
){
    p_adjust_method <- match.arg(p_adjust_method, stats::p.adjust.methods)

    df0 <- stats::model.frame(formula = formula, data = data, drop.unused.levels = TRUE)
    colnames(df0) <- c("y", "x")
    df0 <- df0[ !is.na(df0$y), ]

    n_fct_lvl <- length(unique(df0$x))
    if (n_fct_lvl < 3) stop("Factor levels should be more than 2.")

    if ( ! is_normality(df0, y ~ x) )  # from "./utils.R"
        warning("Data is not normal distribution.")

    if ( levene_test(df0, y ~ x)[["is_var_equal"]] ) # from "./homoscedasticity.R"
        warning("Variations between groups are similar. Please consider Tukey test.")

    # Pre-hoc ====
    pre_hoc <- stats::oneway.test(y ~ x, df0, var.equal = FALSE)
    pre_hoc_pass <- pre_hoc$p.value < alpha

    # Descriptive ====
    desc_mat <- summarize(df0, y ~ x)  # from "./utils.R"
    desc_mat <- desc_mat[order(desc_mat[, "mean"], decreasing = TRUE), ]

    group_means <- desc_mat[, "mean"]
    group_sizes <- desc_mat[, "length"]
    group_names <- names(group_means)
    group_vars <- desc_mat[, "sd"] ^ 2

    is_var_equal <- levene_test(df0, y ~ x)[["is_var_equal"]]

    if (is_var_equal){
        if (length(unique(group_sizes)) == 1)
            warning("Homoscedastic and balanced data, please consider Tukey-HSD test.")
        if (length(unique(group_sizes)) > 1)
            warning("Homoscedastic and unbalanced data, please consider Tukey-Kramer test.")
    }

    # Labels comparison ====
    group_comparisons <- outer2(group_names, function(x1, x2) paste(x1, x2, sep = " |vs| "))

    # Differences of each comparison ====
    group_diff <- outer2(group_means, "-")

    # Standard error ====
    .calc_group_SE <- function(x1, x2){
        var1 <- group_vars[x1]
        var2 <- group_vars[x2]
        n1 <- group_sizes[x1]
        n2 <- group_sizes[x2]
        group_SE <- sqrt( (var1 / n1 + var2 / n2) / 2 )
        return(group_SE)
    }
    group_SE <- outer2(group_names, .calc_group_SE)  # from "utils.R"

    # Pooled degree of freedom ====
    .calc_pooled_df <- function(x1, x2)
    {
        var1 <- group_vars[x1]
        var2 <- group_vars[x2]
        n1 <- group_sizes[x1]
        n2 <- group_sizes[x2]
        numerator <- ( (var1 / n1) + (var2 / n2) ) ^ 2
        denominator <- ( ((var1 / n1) ^ 2) / (n1 - 1) ) + ( ((var2 / n2) ^ 2) / (n2 - 1) )
        pooled_df <- numerator / denominator
        return(pooled_df)
    }
    group_pooled_df <- outer2(group_names, .calc_pooled_df)  # from "utils.R"

    # q-values (Studendized range) ====
    group_qvals <- abs(group_diff / group_SE)

    # p-values ====
    group_pvals <- stats::ptukey(
        q = group_qvals,
        nmeans = length(group_names),
        df = group_pooled_df,
        lower.tail = FALSE
    )
    group_padjs <- stats::p.adjust(group_pvals, method = p_adjust_method)

    # Critical q-values ====
    q_crit <- stats::qtukey(
        p = alpha,
        nmeans = length(group_names),
        df = group_pooled_df,
        lower.tail = FALSE
    )

    # Confidence interval ====
    group_confint <- q_crit * group_SE
    diff_lwr <- group_diff - group_confint
    diff_upr <- group_diff + group_confint

    # Compact letter display ====
    group_cld <- compact_letter_display(
        groups = group_names,
        means = group_means,
        comparisons = group_comparisons,
        pval = group_padjs,
        alpha = alpha,
        descending = TRUE
    )

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Output ====
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    desc_df <- data.frame(
        row.names = NULL,
        GROUP = group_names,
        N = group_sizes,
        AVG = group_means,
        SD = desc_mat[, "sd"],
        MED = desc_mat[, "median"],
        MIN = desc_mat[, "min"],
        MAX = desc_mat[, "max"],
        CLD = group_cld
    )

    comparisons_df <- data.frame(
        comparisons = group_comparisons,
        diff = group_diff,
        lwr = diff_lwr,
        upr = diff_upr,
        StdErr = group_SE,
        qval = group_qvals,
        qcrit = q_crit,
        pval = group_pvals,
        padj = group_padjs
    )

    ret <- list(
        tests = "Welch's ANOVA + Games-Howell test",
        pre_hoc = pre_hoc,
        post_hoc = comparisons_df,
        cld = desc_df
    )

    return(ret)

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ## Testing ====
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    set.seed(1)
    df0 <- data.frame(
        group = rep(c("A", "B", "C"), each = 10),
        norm_data = c(rnorm(10, -1, 3), rnorm(10, 3, 1), rnorm(10, 0, 1)),
        skew_data = c(sqrt(rnorm(10, -7.2, 2) ^ 2), runif(10), runif(10))
    )

    out <- Games_Howell_test(df0, norm_data ~ group)
    out$cld
}






