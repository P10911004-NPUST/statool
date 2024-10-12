# Fisher Least Significant Difference (Fisher-LSD) test


Fisher_LSD_test <- function(
        data,
        formula,
        alpha = 0.05,
        p_adjust_method = "none",
        descending = TRUE
){
    if (!isTRUE(descending) & !isFALSE(descending)) descending <- TRUE
    p_adjust_method <- match.arg(p_adjust_method, stats::p.adjust.methods)

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

    desc_mat <- desc_mat[order(desc_mat[, "mean"], decreasing = descending), ]

    group_means <- desc_mat[, "mean"]
    group_sizes <- desc_mat[, "length"]
    group_names <- names(group_means)


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
    if (length(unique(group_sizes)) > 1) {
        # Unbalanced data
        .calc_group_SE <- function(x1, x2){
            n1 <- group_sizes[x1]
            n2 <- group_sizes[x2]
            res <- sqrt( MSE * ( (1 / n1) + (1 / n2) ) )
            return(res)
        }
        group_SE <- outer2(group_names, .calc_group_SE)
    } else {
        # Balanced data
        group_SE <- sqrt( MSE / mean(group_sizes) )
    }

    # t-values (Studendized range) ====
    group_tvals <- abs(group_diff / group_SE)

    # p-values ====
    group_pvals <- stats::pt(
        q = group_tvals,
        df = DFerror,
        lower.tail = FALSE
    )
    group_padjs <- stats::p.adjust(group_pvals, method = p_adjust_method)

    # Critical t-values ====
    t_crit <- stats::qt(
        p = alpha / 2,
        df = DFerror,
        lower.tail = FALSE
    )

    # Confidence interval ====
    group_confint <- t_crit * group_SE
    diff_lwr <- group_diff - group_confint
    diff_upr <- group_diff + group_confint
    avg_lwr <- group_means - group_confint
    avg_upr <- group_means + group_confint

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
        LCL = avg_lwr,
        UCL = avg_upr,
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
        tval = group_tvals,
        tcrit = t_crit,
        pval = group_pvals,
        padj = group_padjs
    )


    res <- list(
        result = desc_df,
        comparisons = comparisons_df
    )

    return(res)
}





