Dunnett_test <- function(
        data,
        formula,
        control,
        alpha           = 0.05,
        p_adjust_method = "none",
        nsim            = 1e5,
        seed            = 1,
        descending      = NULL
){
    set.seed(seed)

    df0 <- stats::model.frame(formula, data)
    colnames(df0) <- c("y", "x")
    df0 <- df0[!is.na(df0$y), ]

    desc_mat <- summarize(df0, y ~ x)
    if ( ! is.null(descending) )
        desc_mat <- desc_mat[order(desc_mat[, "mean"], decreasing = descending), ]

    group_sizes <- desc_mat[, "length"]
    group_means <- desc_mat[, "mean"]
    group_names <- names(group_means)
    treat_names <- setdiff(group_names, control)

    # ANOVA model ====
    aov_mod <- stats::aov(y ~ x, df0)
    pre_hoc <- anova(aov_mod)

    # Degree of freedom within group ====
    DFerror <- stats::df.residual(aov_mod)

    # Mean square error ====
    MSE <- sum(stats::residuals(aov_mod) ^ 2) / DFerror

    group_comparisons <- vapply(
        X = treat_names,
        FUN = function(trt) paste(control, trt, sep = " |vs| "),
        FUN.VALUE = character(1)
    )

    group_diff <- vapply(
        X = treat_names,
        FUN = function(trt) group_means[control] - group_means[trt],
        FUN.VALUE = numeric(1)
    )

    group_SE <- vapply(
        X = treat_names,
        FUN = function(trt) sqrt((MSE / group_sizes[control]) + (MSE / group_sizes[trt])),
        FUN.VALUE = numeric(1)
    )

    # t-statistics for each treatment groups ====
    group_tvals <- vapply(
        X = treat_names,
        FUN = function(trt) {
            n0 <- desc_mat[rownames(desc_mat) == control, "length"]
            n1 <- desc_mat[rownames(desc_mat) == trt, "length"]
            avg0 <- desc_mat[rownames(desc_mat) == control, "mean"]
            avg1 <- desc_mat[rownames(desc_mat) == trt, "mean"]
            abs( (avg0 - avg1) / sqrt(MSE * ((1 / n0) + (1 / n1))) )
        },
        FUN.VALUE = numeric(1)
    )

    # Calculate correlation matrix
    .calculate_corr_mat <- function(i, j){
        n0 <- 1 / group_sizes[control]
        ni <- 1 / group_sizes[i]
        nj <- 1 / group_sizes[j]
        n0 / sqrt((ni + n0) * (nj + n0))
    }
    corr_mat <- outer(treat_names, treat_names, .calculate_corr_mat)
    diag(corr_mat) <- 1

    # Simulate multivariate t-distribution
    L <- chol(corr_mat)
    Z <- matrix(stats::rnorm(length(treat_names) * nsim), nsim) %*% L
    s <- sqrt(stats::rchisq(nsim, DFerror) / DFerror)
    t_sim <- Z / s

    # Calculate maximum absolute values
    max_abs_t <- apply(abs(t_sim), 1, max)

    # Calculate critical value and p-values
    t_crit <- stats::quantile(max_abs_t, 1 - alpha, names = FALSE)
    pvals <- vapply(
        X = group_tvals,
        FUN = function(t) mean(max_abs_t >= t),
        FUN.VALUE = numeric(1)
    )
    padjs <- stats::p.adjust(pvals, method = p_adjust_method)

    group_confint <- t_crit * group_SE

    df_comparisons <- data.frame(
        comparisons = group_comparisons,
        diff = group_diff,
        lwr = group_diff - group_confint,
        upr = group_diff + group_confint,
        StdErr = group_SE,
        tval = group_tvals,
        tcrit = t_crit,
        pval = pvals,
        padj = padjs,
        "signif." = pval2asterisk(padjs)
    )

    asterisks <- stats::setNames(
        object = c("", df_comparisons[["signif."]]),
        nm = c(control, rownames(df_comparisons))
    )
    asterisks <- asterisks[match(names(asterisks), group_names)]
    rownames(df_comparisons) <- NULL

    desc_df <- data.frame(
        row.names = NULL,
        GROUP = group_names,
        N = group_sizes,
        AVG = group_means,
        SD = desc_mat[, "sd"],
        MED = desc_mat[, "median"],
        MIN = desc_mat[, "min"],
        MAX = desc_mat[, "max"],
        CLD = asterisks
    )

    ret <- list(
        tests = "Dunnett's test",
        pre_hoc = pre_hoc,
        post_hoc = df_comparisons,
        cld = desc_df
    )

    return(ret)

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ## Testing ====
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    set.seed(1)
    df0 <- data.frame(
        fertilizer = rep(c("A", "B", "C", "D", "E"), each = 10),
        yield = c(
            stats::rnorm(10, 8, 1.2),
            stats::rnorm(10, 9.1, 1.1),
            stats::rnorm(10, 8.5, 1.3),
            stats::rnorm(10, 10, 1.05),
            stats::rnorm(10, 7, 1.17))
    )

    control <- subset(df0, fertilizer == "A")$y
    treatment <- subset(df0, fertilizer != "A")
    treatment <- dataframe_to_list(
        data = treatment,
        formula = yield ~ fertilizer
    )

    out <- Dunnett_test(df0, yield ~ fertilizer, control = "A")
    out$post_hoc
    DescTools::DunnettTest(yield ~ fertilizer, df0, control = "A")

}


