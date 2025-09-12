#' @title
#' Student's t-Test
#'
#' @description
#' Performs one and two sample t-tests on vectors of data.
#'
#'
#' @param data A data.frame
#' @param formula A formula y ~ x specifying the response and independant variable.
#' @param y1 A numeric vector of data values
#' @param y2 A numeric vectoe of another data values
#' @param alternative a character string specifying the alternative hypothesis,
#'      must be one of "two.sided" (default), "greater" or "less".
#' @param var_equal A logical variable indicating whether to treat the two variances as being equal.
#'      If TRUE then the pooled variance is used to estimate the variance otherwise the Welch (or Satterthwaite) approximation to the degrees of freedom is used.
#' @param alpha Confidence level (1 - alpha) of the interval.
#' @param paired A logical indicating whether you want a paired t-test.
#'
#' @returns A data.frame
#' @export
#'
#' @examples
#' y1 <- rnorm(10)
#' y2 <- rnorm(10, 5, 3)
#' t_test(y1 = y1, y2 = y2)
t_test <- function(
        data,
        formula = NULL,
        y1 = NULL,
        y2 = NULL,
        alternative = c("two.sided", "less", "greater"),
        var_equal = NULL,
        alpha = 0.05,
        paired = FALSE
) {
    if (is.null(formula))
    {
        if (inherits(y1, "character"))
            y1 <- data[[y1]]

        if (inherits(y2, "character"))
            y2 <- data[[y2]]

        if (is.null(var_equal))
        {
            fstats <- F_test(data, x = y1, y = y2)
            var_equal <- fstats[["pval"]] > 0.05
        }

        tstats <- stats::t.test(
            x = y1,
            y = y2,
            var.equal = var_equal,
            alternative = alternative,
            conf.level = 1 - alpha
        )

    } else {

        if (is.null(var_equal))
        {
            fstats <- F_test(data, formula)
            var_equal <- fstats[["pval"]] > 0.05
        }

        tstats <- stats::t.test(
            formula,
            data = data,
            var.equal = var_equal,
            alternative = alternative,
            conf.level = 1 - alpha
        )
    }

    tstats <- data.frame(
        comparison = tstats[["data.name"]],
        t = unname(tstats[["statistic"]]),
        df = unname(tstats[["parameter"]]),
        pval = tstats[["p.value"]],
        StdErr = tstats[["stderr"]],
        conf_int_lower = tstats[["conf.int"]][1],
        conf_int_upper = tstats[["conf.int"]][2]
    )

    return(tstats)
}
