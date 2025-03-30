#' Dunn's test
#'
#' Represent significance statements resulting from all-pairwise comparisons.
#'
#' @param data A data frame in which the variables specified in the formula will be found.
#' @param formula A formula specifying the model.
#' @param alpha Numeric value range from 0 to 1 (Default: 0.05).
#' @param p_adjust_method Character string (default: "holm"). Other options: `stats::p.adjust.methods`.
#'
#' @return A list with two vectors. \
#' 1. result: consists of descriptive statistics and compact letter display;
#' 2. comparisons: includes statistics parameters for each pairwise comparisons.
#'
#' @export
#' @author Joon-Keat Lai
#'
#' @references
#' Dunn, O.J., 1964.
#' Multiple Comparisons Using Rank Sums.
#' Technometrics 6, 241–252.
#' https://doi.org/10.1080/00401706.1964.10490181
#'
#' @examples
#' set.seed(1)
#' df0 <- data.frame(
#'     group = rep(c("A", "B", "C"), each = 10),
#'     norm_data = c(rnorm(10, -1, 2), rnorm(10, 3, 2), rnorm(10, 0, 1.5)),
#'     skew_data = c(sqrt(rnorm(10, -7.2, 2) ^ 2), runif(10), runif(10))
#' )
#'
#' out <- Dunn_test(df0, skew_data ~ group)
#' res <- out$result
#' res
#' #>   GROUPS  N  AVG       SD  MED MIN MAX CLD
#' #> 1      A 10 25.5 3.027650 25.5  21  30   a
#' #> 2      C 10 12.5 6.258328 12.5   1  20   b
#' #> 3      B 10  8.5 5.082650  7.5   2  16   b
Dunn_test <- function(
        data,
        formula,
        alpha = 0.05,
        p_adjust_method = "holm"
){
    p_adjust_method <- match.arg(p_adjust_method, stats::p.adjust.methods)
    adjusted_alpha <- 1 - stats::p.adjust(1 - alpha, p_adjust_method)

    df0 <- stats::model.frame(formula = formula, data = data, drop.unused.levels = TRUE)
    colnames(df0) <- c("y", "x")
    df0 <- df0[ !is.na(df0$y), ]

    if (is_normality(df0, y ~ x))  # is_normality() <<< ./normality.R
        warning("The response variable follows a normal distribution.")

    n_fct <- length(unique(df0$x))
    if (n_fct < 3) stop("Factor levels should be more than 2.")

    df0$y <- rank(df0$y)

    desc_mat <- with(
        data = df0,
        expr = vapply(
            X = c("sum", "length", "mean", "sd", "median", "min", "max"),
            FUN = function(fns) tapply(y, x, fns),
            FUN.VALUE = numeric(n_fct)
        )
    )
    desc_mat <- desc_mat[order(desc_mat[, "mean"], decreasing = TRUE), ]

    group_means <- desc_mat[, "mean"]
    group_sizes <- desc_mat[, "length"]
    group_names <- names(group_means)
    group_comparisons <- outer2(group_names, function(x1, x2) paste(x1, x2, sep = " |vs| "))
    group_n <- outer2(group_sizes, function(x1, x2) paste(x1, x2, sep = "|"))
    group_diffs <- stats::setNames(outer2(group_means, "-"), group_comparisons) # outer2() <<< ./utils.R
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

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ## Compact letter display ====
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    group_cld <- compact_letter_display(
        # compact_letter_display() <<< ./compact_letter_display.R
        groups = group_names,
        means = group_means,
        comparisons = group_comparisons,
        pval = group_padjs,
        alpha = alpha,
        descending = TRUE
    )

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ## Output ====
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
        row.names = NULL,
        comparisons = group_comparisons,
        n = group_n,
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

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ## Testing ====
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    set.seed(1)
    df0 <- data.frame(
        group = rep(c("A", "B", "C"), each = 10),
        norm_data = c(rnorm(10, -1, 2), rnorm(10, 3, 2), rnorm(10, 0, 1.5)),
        skew_data = c(sqrt(rnorm(10, -7.2, 2) ^ 2), runif(10), runif(10))
    )

    out <- Dunn_test(df0, norm_data ~ group)
    res <- out$result
    res

}





