#' Tukey Honestly Significant Difference (Tukey-HSD) test
#'
#' Represent significance statements resulting from all-pairwise comparisons.
#'
#' @param data A data frame in which the variables specified in the formula will be found.
#' @param formula A formula specifying the model.
#' @param alpha Numeric value range from 0 to 1 (default: 0.05). The error tolerance.
#' @param p_adjust_method
#' A character string indicating which method to use for p-value adjustment (default: "none").
#' For Tukey-HSD, this option should always be "none", and the calculated padj is "Tukey-adjusted p-value".
#' Other options: `stats::p.adjust.methods` (NOT recommended).
#'
#' @return A list with two vectors.
#' 1. result: consists of descriptive statistics and compact letter display;
#' 2. comparisons: includes statistics parameters for each pairwise comparisons.
#'
#' @export
#' @author Joon-Keat Lai
#'
#' @references
#' Howell, D.C. 2013.
#' Statistical methods for psychology (8th ed.). pg. 393.
#' Wadsworth Cengage Learning, Belmont, CA.
#'
#' Tukey, John W. 1949.
#' Comparing Individual Means in the Analysis of Variance.
#' Biometrics, vol. 5, no. 2, pp. 99-114. JSTOR.
#' https://doi.org/10.2307/3001913.
#'
#' @examples
#' set.seed(1)
#' df0 <- data.frame(
#'     group = as.factor(rep(c("A", "B", "C"), each = 10)),
#'     norm_data = c(rnorm(10, -1, 2), rnorm(10, 3, 2), rnorm(10, 0, 2)),
#'     skew_data = c(sqrt(rnorm(10, -7.2, 2) ^ 2), runif(10), runif(10))
#' )
#'
#' out <- Tukey_HSD_test(df0, norm_data ~ group)
#' out$cld
#' #>   GROUP  N        AVG       SD         MED       MIN      MAX CLD
#' #> 1     B 10  3.4976899 2.139030  3.98374456 -1.429400 6.023562   a
#' #> 2     C 10 -0.2673465 1.911215  0.01843624 -3.978703 1.837955   b
#' #> 3     A 10 -0.7355944 1.561172 -0.48684890 -2.671257 2.190562   b
Tukey_HSD_test <- function(
        data,
        formula,
        alpha = 0.05,
        p_adjust_method = "none"
){
    p_adjust_method <- match.arg(p_adjust_method, choices = stats::p.adjust.methods)

    df0 <- stats::model.frame(formula = formula, data = data, drop.unused.levels = TRUE)
    colnames(df0) <- c("y", "x")
    df0 <- df0[ !is.na(df0$y), ]

    n_fct_lvl <- length(unique(df0$x))
    if (n_fct_lvl < 3) warning("Factor levels should be more than 2.")

    # is_normality(), is_unbalance() <<< ./utils.R
    if ( ! is_normality(df0, y ~ x) ) warning("Data is not normal distribution.")
    if ( is_unbalance(df0, y ~ x) ) warning("Unbalanced design. Please consider Tukey-Kramer test")
    # levene_test() <<< ./homoscedasticity.R
    if ( ! levene_test(df0, y ~ x)[["is_var_equal"]] )
        warning("Variations between groups are different. Please consider Games-Howell test.")

    # Pre-hoc ====
    pre_hoc <- stats::oneway.test(y ~ x, df0, var.equal = TRUE)
    pre_hoc_pass <- pre_hoc$p.value < alpha

    # Descriptive ====
    desc_mat <- with(
        data = df0,
        expr = vapply(
            X = c("length", "mean", "sd", "median", "min", "max"),
            FUN = function(fns) tapply(y, x, fns),
            FUN.VALUE = numeric(n_fct_lvl)
        )
    )

    desc_mat <- desc_mat[order(desc_mat[, "mean"], decreasing = TRUE), ]

    group_means <- desc_mat[, "mean"]
    group_sizes <- desc_mat[, "length"]
    group_names <- names(group_means)

    if (length(unique(group_sizes)) > 1)
        warning("Unbalanced data, please consider Tukey-Kramer test.")

    # ANOVA model ====
    aov_mod <- stats::aov(formula = y ~ x, data = df0)

    # Degree of freedom (within group) ====
    DFerror <- stats::df.residual(aov_mod)

    # Mean square error ====
    MSE <- sum(stats::residuals(aov_mod) ^ 2) / DFerror

    # Labels comparison ====
    group_comparisons <- outer2(group_names, function(x1, x2) paste(x1, x2, sep = " |vs| "))

    # Differences of each comparison ====
    group_diff <- outer2(group_means, "-")

    # Standard error ====
    group_SE <- sqrt( MSE / mean(group_sizes) )

    # q-values (Studendized range) ====
    group_qvals <- abs(group_diff / group_SE)

    # p-values ====
    group_pvals <- stats::ptukey(
        q = group_qvals,
        nmeans = length(group_names),
        df = DFerror,
        lower.tail = FALSE
    )
    group_padjs <- stats::p.adjust(group_pvals, method = p_adjust_method)

    # Critical q-values ====
    q_crit <- stats::qtukey(
        p = alpha,
        nmeans = length(group_names),
        df = DFerror,
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

    # Output ====
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
        tests = "Fisher's ANOVA + Tukey-HSD test",
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
        group = as.factor(rep(c("A", "B", "C"), each = 10)),
        norm_data = c(rnorm(10, -1, 2), rnorm(10, 3, 2), rnorm(10, 0, 2)),
        skew_data = c(sqrt(rnorm(10, -7.2, 2) ^ 2), runif(10), runif(10))
    )

    out <- Tukey_HSD_test(df0, norm_data ~ group)
    out$cld
    # stats::TukeyHSD(aov(norm_data ~ group, df0))
}




