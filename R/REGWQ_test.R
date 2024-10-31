#' Ryan, Einot, Gabriel, Welsh Studentized Range Q (REGWQ) test
#'
#' Represent significance statements resulting from all-pairwise comparisons.
#'
#' @param data A data frame in which the variables specified in the formula will be found.
#' @param formula A formula specifying the model.
#' @param alpha Numeric value range from 0 to 1 (default: 0.05). The error tolerance.
#' @param p_adjust_method A character string (default: "none"). Options please refer to `stats::p.adjust.methods`.
#'
#' @return A list with two vectors. \
#' 1. result: consists of descriptive statistics and compact letter display;
#' 2. comparisons: includes statistics parameters for each pairwise comparisons.
#' @export
#'
#' @examples
#' rawdata <- data.frame(
#'     m1 = c(51, 84, 50, 48, 79, 61, 53, 54),
#'     m2 = c(82, 91, 92, 80, 52, 85, 73, 74),
#'     m3 = c(79, 84, 74, 98, 63, 83, 85, 58),
#'     m4 = c(85, 80, 65, 71, 67, 51, 63, 93),
#'     m5 = c(37, 40, 61, 51, 76, 55, 60, 70)
#' )
#'
#' rawdata <- stats::reshape(
#'     data = rawdata,
#'     direction = "long",
#'     varying = colnames(rawdata),
#'     v.names = "val",
#'     times = colnames(rawdata),
#'     timevar = "group"
#' )
#'
#' stats_res <- REGWQ_test(rawdata, val ~ group)
#' stats_res[["result"]]["CLD"]
#' #>    GROUPS N    AVG       SD  MED MIN MAX CLD
#' #> m2     m2 8 78.625 12.80555 81.0  52  92   a
#' #> m3     m3 8 78.000 12.82854 81.0  58  98   a
#' #> m4     m4 8 71.875 13.47418 69.0  51  93  ab
#' #> m1     m1 8 60.000 13.87701 53.5  48  84   b
#' #> m5     m5 8 56.250 13.51983 57.5  37  76   b
#'
#' @references
#' Howell, D.C. 2013.
#' Statistical methods for psychology (8th ed.). pg. 393.
#' Wadsworth Cengage Learning, Belmont, CA.
REGWQ_test <- function(
        data,
        formula,
        alpha = 0.05,
        p_adjust_method = "none"
){

    p_adjust_method <- match.arg(p_adjust_method, choices = stats::p.adjust.methods)

    y_name <- as.character(formula)[2]
    x_name <- as.character(formula)[3]

    df0 <- data.frame(
        y = data[[y_name]],
        x = data[[x_name]]
    )

    df0 <- df0[!is.na(df0$y), ]

    # Descriptive ====
    desc_mat <- with(
        data = df0,
        expr = sapply(
            X = c("length", "mean", "sd", "median", "min", "max"),
            FUN = function(fns) tapply(y, x, fns)
        )
    )

    desc_mat <- desc_mat[order(desc_mat[, "mean"], decreasing = TRUE), ]

    group_means <- desc_mat[, "mean"]
    group_sizes <- desc_mat[, "length"]
    group_names <- names(group_means)

    group_comparisons <- outer2(group_names, function(x1, x2) paste(x1, x2, sep = " |vs| "))
    group_diff <- outer2(group_means, "-")

    if (length(unique(group_sizes)) > 1)
        stop("Unbalanced data, please consider Tukey-Kramer test.")

    # ANOVA model ====
    aov_mod <- stats::aov(formula = y ~ x, data = df0)

    # Degree of freedom (within group) ====
    DFerror <- stats::df.residual(aov_mod)

    # Mean square error ====
    # same as sum(stats::residuals(aov_mod) ^ 2) / DFerror
    MSE <- stats::deviance(aov_mod) / DFerror

    k <- length(group_names)

    group_SE <- sqrt( MSE / mean(group_sizes) )

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # .calc_pvals() ====
    # Create object within the REGWQ_test() function,
    # otherwise, the super assignment will assign to the global environment
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    out_modified_alpha <- c()
    out_k_subset <- c()
    out_group_qvals <- c()
    out_group_q_crits <- c()

    .calc_pvals <- function(x1, x2){
        x1_mean <- group_means[x1]
        x2_mean <- group_means[x2]

        x1_index <- sapply(x1, function(x) grep(x, group_names))
        x2_index <- sapply(x2, function(x) grep(x, group_names))

        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        ## k_subset ====
        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        k_subset <- x2_index - x1_index + 1

        ### Output k_subset matrix to the outside environment of this function
        if (TRUE){
            out_k_subset <- matrix(k_subset, nrow = k)
            rownames(out_k_subset) <- colnames(out_k_subset) <- group_names
            out_k_subset <<- out_k_subset[upper.tri(out_k_subset)]
        }


        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        ## Modified alpha ====
        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        .modify_alpha <- function(k_subset) {
            ifelse(
                test = ( k_subset < (k - 1) ),
                yes = return(1 - ( (1 - alpha) ^ (k_subset / k) )),
                no = return(alpha)
            )
        }

        modified_alpha <- sapply(k_subset, .modify_alpha)

        ### Output alpha matrix to the outside environment of this function
        if (TRUE) {
            out_modified_alpha <-  matrix(modified_alpha, nrow = k)
            rownames(out_modified_alpha) <- colnames(out_modified_alpha) <- group_names
            out_modified_alpha <<- out_modified_alpha[upper.tri(out_modified_alpha)]
        }


        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        ## Q-vals ====
        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        group_qvals <- abs(x1_mean - x2_mean) / sqrt( MSE / mean(group_sizes) )

        ### Output alpha matrix to the outside environment of this function
        if (TRUE){
            out_group_qvals <-  matrix(group_qvals, nrow = k)
            rownames(out_group_qvals) <- colnames(out_group_qvals) <- group_names
            out_group_qvals <<- out_group_qvals[upper.tri(out_group_qvals)]
        }


        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        ## Q-critical ====
        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        group_q_crits <- suppressWarnings(
            stats::qtukey(
                p = modified_alpha,
                nmeans = k_subset,
                df = DFerror,
                lower.tail = FALSE
            )
        )

        ### Output q_crit matrix to the outside environment of this function
        if (TRUE){
            out_group_q_crits <-  matrix(group_q_crits, nrow = k)
            rownames(out_group_q_crits) <- colnames(out_group_q_crits) <- group_names
            out_group_q_crits <<- out_group_q_crits[upper.tri(out_group_q_crits)]
        }


        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        ## p-value ====
        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        group_pvals <- suppressWarnings(
            stats::ptukey(
                q = group_qvals,
                nmeans = k_subset,
                df = DFerror,
                lower.tail = FALSE
            )
        )

        return(group_pvals)
    }  # End of .calc_pvals()


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Iterate the p-values ====
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    group_pvals <- outer(group_names, group_names, .calc_pvals)
    group_pvals[is.nan(group_pvals)] <- 1

    # Change the p-values row-wise and elememt-wise
    tmp_mat <- group_pvals
    n_rows <- nrow(group_pvals)
    n_cols <- ncol(group_pvals)
    for (row in 1:n_rows){
        for (col in 1:(n_cols - 1)){
            tmp_mat[row, col] <- ifelse(
                test = (group_pvals[row, col + 1] > alpha) & (group_pvals[row, col] <= alpha),
                yes = 1,
                no = group_pvals[row, col]
            )
        }
    }
    group_pvals <- tmp_mat
    group_pvals <- group_pvals[upper.tri(group_pvals)]

    group_padjs <- stats::p.adjust(group_pvals, method = p_adjust_method)

    # Compact letter display ====
    group_cld <- compact_letter_display(
        groups = group_names,
        means = group_means,
        comparisons = group_comparisons,
        pval = group_padjs,
        alpha = alpha,
        descending = TRUE
    )

    # Confidence interval ====
    group_confint <- out_group_q_crits * group_SE
    diff_lwr <- group_diff - group_confint
    diff_upr <- group_diff + group_confint


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Output ====
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
        k_subset = out_k_subset,
        alpha_k = out_modified_alpha,
        diff = group_diff,
        diff_lwr = diff_lwr,
        diff_upr = diff_upr,
        StdErr = group_SE,
        qval = out_group_qvals,
        qcrit = out_group_q_crits,
        pval = group_pvals,
        padj = group_padjs
    )


    res <- list(
        result = desc_df,
        comparisons = comparisons_df
    )

    return(res)
}








