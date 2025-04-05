#' One-way aligned rank transform ANOVA
#'
#' @param data A data frame.
#' @param formula A formula specifying the model, for example, yield ~ fertilizer_level.
#' @param alpha Numeric value range from 0 to 1 (default: 0.05). The error tolerance.
#' @param p_adjust_method A character string (default: "none"). Other options: `stats::p.adjust.methods`.
#'
#' @importFrom car Anova
#'
#' @return ANOVA (type 3) result
#' @export
#'
#' @author Joon-Keat Lai
#'
#' @references
#' Wobbrock, J.O., Findlater, L., Gergle, D. and Higgins, J.J. (2011).
#' The aligned rank transform for nonparametric factorial analyses using only ANOVA procedures.
#' Proceedings of the ACM Conference on Human Factors in Computing Systems (CHI '11).
#' Vancouver, British Columbia (May 7-12, 2011). New York: ACM Press, pp. 143-146.
#' https://doi.org/10.1145/1978942.1978963
#'
#' Elkin, L.A., Kay, M., Higgins, J.J. and Wobbrock, J.O. (2021).
#' An aligned rank transform procedure for multifactor contrast tests.
#' Proceedings of the ACM Symposium on User Interface Software and Technology (UIST '21).
#' Virtual Event (October 10-14, 2021). New York: ACM Press, pp. 754-768.
#' https://doi.org/10.1145/3472749.3474784
#'
#' @examples
#' set.seed(1)
#' df0 <- data.frame(
#'     group = as.factor(rep(c("A", "B", "C"), each = 10)),
#'     skew_data1 = c(rlnorm(10, -1, 2), sqrt(rcauchy(10, 3, 2) ^ 2), rlnorm(10, 0, 1.5)),
#'     skew_data2 = c(sqrt(rnorm(10, -7.2, 2) ^ 2), runif(10), runif(10))
#' )
#' out <- art_1(df0, skew_data1 ~ group)
#' print(out$pre_hoc)
#' #> Anova Table (Type 3)
#' #> art(skew_data1) ~ group
#' #>             Sum Sq Df F value    Pr(>F)
#' #> (Intercept)  883.6  1  16.002 0.0004423 ***
#' #> x            756.6  2   6.851 0.0039229 **
#' #> Residuals   1490.9 27
#' #> ---
#' #> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#'
#' print(out$cld)
#' #>   GROUP  N      AVG       SD       MED        MIN      MAX CLD
#' #> 1     B 10 4.547070 2.543874 4.2027752 0.87471697 8.349396   a
#' #> 2     C 10 2.279905 1.429011 2.4855260 0.05058867 4.119593  ab
#' #> 3     A 10 1.437723 2.687525 0.6211093 0.06916521 8.940233   b

art_1 <- function(data, formula, alpha = 0.05, p_adjust_method = "none")
{
    df0 <- stats::model.frame(formula, data)
    y_name <- colnames(df0)[1]
    x_name <- colnames(df0)[2]

    colnames(df0) <- c("y", "x")
    df0 <- df0[!is.na(df0$y), ]

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ## Pre-hoc ====
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    aov_mod <- stats::aov(y ~ x, df0)
    df0$residual <- stats::residuals(aov_mod)

    group_means <- with(df0, tapply(y, x, mean))
    estimated_effect <- group_means - mean(df0$y)

    df0$estimated_effect <- vapply(
        X = df0$x,
        FUN = function(x) estimated_effect[match(x, names(estimated_effect))],
        FUN.VALUE = numeric(1)
    )

    df0$aligned_Y <- with(df0, residual + estimated_effect)  # aligned response
    df0$ranked_Y <- rank(df0$aligned_Y, ties.method = "average")  # ranked response

    aov_mod <- stats::aov(ranked_Y ~ x, df0)
    pre_hoc <- car::Anova(aov_mod, type = 3)
    attr(pre_hoc, "heading") <- sprintf("Anova Table (Type 3)\nart(%s) ~ %s", y_name, x_name)

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ## Post-hoc ====
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if (!is_normality(df0, ranked_Y ~ x))  # is_normality() <<< normality.R
        warning("Aligned rank transform may be appropriate for this dataset")

    is_balanced <- length(unique(table(df0$x))) == 1
    if (is_balanced) {
        tests <- "ART + Tukey-HSD on ART-response"
        tukey <- Tukey_HSD_test(df0, ranked_Y ~ x, alpha)
    } else {
        tests <- "ART + Tukey-Kramer on ART-response"
        tukey <- Tukey_Kramer_test(df0, ranked_Y ~ x, alpha)
    }

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ## DO NOT use ranked data in the descriptive table (use `y` rather than `ranked_y`)
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    desc_mat <- summarize(df0, y ~ x)
    desc_mat <- desc_mat[match(tukey$cld$GROUP, rownames(desc_mat)), ]

    tukey$cld <- data.frame(
        row.names = NULL,
        GROUP = tukey$cld$GROUP,
        N = tukey$cld$N,
        AVG = desc_mat[, "mean"],
        SD = desc_mat[, "sd"],
        MED = desc_mat[, "median"],
        MIN = desc_mat[, "min"],
        MAX = desc_mat[, "max"],
        CLD = tukey$cld$CLD
    )

    colnames(df0)[colnames(df0) == "y"] <- y_name
    colnames(df0)[colnames(df0) == "aligned_Y"] <- "Y'"
    colnames(df0)[colnames(df0) == "ranked_Y"] <- "Y''"
    colnames(df0)[colnames(df0) == "x"] <- x_name

    ret <- list(
        tests = tests,
        data = df0,
        pre_hoc = tukey$pre_hoc,
        post_hoc = tukey$post_hoc,
        cld = tukey$cld
    )

    return(ret)

    ##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ## Test
    ##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    set.seed(1)
    df0 <- data.frame(
        group = as.factor(rep(c("A", "B", "C"), each = 10)),
        skew_data1 = c(rlnorm(10, -1, 2), sqrt(rcauchy(10, 3, 2) ^ 2), rlnorm(10, 0, 1.5)),
        skew_data2 = c(sqrt(rnorm(10, -7.2, 2) ^ 2), runif(10), runif(10))
    )
    out <- art_1(df0, skew_data1 ~ group)
    print(out$pre_hoc)
    print(out$cld)
}


