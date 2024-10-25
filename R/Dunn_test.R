#' Dunn's test
#'
#' Represent significance statements resulting from all-pairwise comparisons.
#'
#' @param data A data frame in which the variables specified in the formula will be found.
#' @param formula A formula specifying the model.
#' @param alpha Numeric value range from 0 to 1 (default: 0.05). The error tolerance.
#' @param p_adjust_method A character string (default: "none"). Options please refer to `stats::p.adjust.methods`.
#' @param descending Logical (default: TRUE). Sort the ranks in descending order.
#'
#' @return A list with two vectors. \
#' 1. result: consists of descriptive statistics and compact letter display;
#' 2. comparisons: includes statistics parameters for each pairwise comparisons.
#' @export

Dunn_test <- function(
        data,
        formula,
        alpha = 0.05,
        p_adjust_method = "holm",
        descending = TRUE
){
    p_adjust_method <- match.arg(p_adjust_method, stats::p.adjust.methods)
    adjusted_alpha <- 1 - stats::p.adjust(1 - alpha, p_adjust_method)

    df0 <- stats::model.frame(formula = formula, data = data, drop.unused.levels = TRUE)
    colnames(df0) <- c("y", "x")
    df0 <- df0[ !is.na(df0$y), ]

    df0$y <- rank(df0$y)

    desc_mat <- with(
        data = df0,
        expr = sapply(
            X = c("sum", "length", "mean", "sd", "median", "min", "max"),
            FUN = function(fns) tapply(y, x, fns)
        )
    )
    desc_mat <- desc_mat[order(desc_mat[, "mean"], decreasing = descending), ]

    group_means <- desc_mat[, "mean"]
    group_sizes <- desc_mat[, "length"]
    group_names <- names(group_means)
    group_comparisons <- outer2(group_names, function(x1, x2) paste(x1, x2, sep = " |vs| "))
    group_diffs <- stats::setNames(outer2(group_means, "-"), group_comparisons)
    n <- sum(group_sizes)

    .ties_correction <- function(x){
        ties_num <- table(x)
        res <- sum( ties_num ^ 3 - ties_num ) / ( 12 * (n - 1) )
        return(res)
    }

    block_A <- ( n * (n + 1) ) / 12
    block_B <- .ties_correction(df0$y)
    block_C <- outer2( group_sizes, function(x1, x2) {( (1 / x1) + (1 / x2) )} )

    group_SEs <- sqrt( (block_A - block_B) * block_C )
    names(group_SEs) <- group_comparisons

    group_zvals <- abs(group_diffs) / group_SEs
    z_crit <- stats::qnorm(adjusted_alpha / 2, lower.tail = FALSE)

    group_pvals <- 2 * stats::pnorm(group_zvals, lower.tail = FALSE)
    group_padjs <- stats::p.adjust(group_pvals, method = p_adjust_method)

    group_diff_thresh <- z_crit * group_SEs

    # Compact letter display ====
    group_cld <- compact_letter_display(
        groups = group_names,
        means = group_means,
        comparisons = group_comparisons,
        pval = group_padjs,
        alpha = alpha,
        descending = descending
    )

    # Output ====
    desc_df <- data.frame(
        GROUPS = group_names,
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
        diff = group_diffs,
        StdErr = group_SEs,
        zval = group_zvals,
        z_crit = z_crit,
        pval = group_pvals,
        padj = group_padjs
    )


    res <- list(
        result = desc_df,
        comparisons = comparisons_df
    )

    return(res)
}


# rawdata <- data.frame(
#     New = c(81, 32, 42, 62, 37, 44, 38, 47, 49, 41),
#     Old = c(48, 31, 25, 22, 30, 30, 32, 15, 40, NA),
#     Control = c(18, 49, 33, 19, 24, 17, 48, 22, NA, NA)
# )
#
# rawdata <- stats::reshape(
#     data = rawdata,
#     direction = "long",
#     varying = colnames(rawdata),
#     v.names = "val",
#     timevar = "group",
#     times = colnames(rawdata),
# )
#
# rawdata <- rawdata[!is.na(rawdata$val), ]
#
# Dunn_test(rawdata, val ~ group)



