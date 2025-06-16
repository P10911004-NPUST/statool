#' Dunn's test
#'
#' Represent significance statements resulting from all-pairwise comparisons.
#'
#' @param data A data frame in which the variables specified in the formula will be found.
#' @param formula A formula specifying the model.
#' @param alpha Numeric value range from 0 to 1 (Default: 0.05).
#' @param p_adjust_method p-value adjustment method (Default: "holm"). See `stats::p.adjust.methods`.
#'
#' @return A list with four vectors.
#' 1. tests: A message showing the statistical methods applied on the dataset.
#' 2. pre_hoc: The result of pmnibus test.
#' 3. post_hoc: includes statistics parameters for each pairwise comparisons.
#' 4. cld: A dataframe reporting the descriptive stats and compact letter display.
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
#' cld <- out$cld
#' cld
#' #>   GROUP  N       AVG        SD       MED        MIN       MAX CLD
#' #> 1     A 10 6.9585396 1.6171293 7.3131184 4.48264090 9.9541191   a
#' #> 2     C 10 0.5674507 0.2713306 0.6236108 0.05893438 0.8762692   b
#' #> 3     B 10 0.4053906 0.2437230 0.3626733 0.12169192 0.7570871   b
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

    if (is_normality(df0, y ~ x))  # from "./normality.R"
        warning("The response variable follows a normal distribution.")

    n_fct_lvl <- length(unique(df0$x))
    if (n_fct_lvl < 3) stop("Factor levels should be more than 2.")

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ## Pre-hoc ====
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    pre_hoc <- with(df0, stats::kruskal.test(y, x))
    pre_hoc_pass <- pre_hoc[["p.value"]] < alpha

    df0$ranked_y <- rank(df0$y)

    desc_mat <- summarize(df0, ranked_y ~ x)  # from "./utils.R"
    desc_mat <- desc_mat[order(desc_mat[, "mean"], decreasing = TRUE), ]

    group_means <- desc_mat[, "mean"]
    group_sizes <- desc_mat[, "length"]
    group_names <- names(group_means)

    group_comparisons <- outer2(   # from "utils.R"
        group_names,
        function(x1, x2) paste(x1, x2, sep = " |vs| ")
    )
    group_n <- outer2(group_sizes, function(x1, x2) paste(x1, x2, sep = "|"))
    group_diffs <- stats::setNames(outer2(group_means, "-"), group_comparisons)

    total_N <- sum(group_sizes)

    .ties_correction <- function(x)
    {
        ties_num <- table(x)
        res <- sum( ties_num ^ 3 - ties_num ) / ( 12 * (total_N - 1) )
        return(res)
    }

    block_A <- ( total_N * (total_N + 1) ) / 12
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
    group_cld <- compact_letter_display(  # from "./compact_letter_display.R"
        groups = group_names,
        means = group_means,
        comparisons = group_comparisons,
        pval = group_padjs,
        alpha = alpha,
        descending = TRUE
    )

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ## DO NOT use ranked data in the descriptive table (use `y` rather than `ranked_y`)
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    desc_mat2 <- summarize(df0, y ~ x)  # from "./utils.R"
    desc_mat2 <- desc_mat2[match(group_names, rownames(desc_mat2)), ]

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ## Output ====
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    desc_df <- data.frame(
        row.names = NULL,
        GROUP = group_names,
        N = group_sizes,
        AVG = desc_mat2[, "mean"],
        SD = desc_mat2[, "sd"],
        MED = desc_mat2[, "median"],
        MIN = desc_mat2[, "min"],
        MAX = desc_mat2[, "max"],
        CLD = group_cld
    )

    if ( ! pre_hoc_pass ) desc_df$CLD <- "a"

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
        tests = "Kruskal-Wallis test + Dunn's test",
        pre_hoc = pre_hoc,
        post_hoc = comparisons_df,
        cld = desc_df
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

    out <- Dunn_test(df0, skew_data ~ group)
    cld <- out$cld
    cld
}





