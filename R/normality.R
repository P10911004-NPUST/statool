# Normality test ====

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Anderson-Darling ====
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#' Anderson-Darling normality test
#'
#' Test the data normality (normal distribution)
#'
#' @param x A numeric vector
#'
#' @return A list with three elements
#' 1. is_normality: logical value. TRUE: x is normal distribution
#' 2. statistic: the Anderson-Darling test statistic value.
#' The higher the value, the lower the probability of normality.
#' 3. pval: the p-value of the statistics. pval > 0.05 indicates normal distribution.
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- stats::rnorm(100)  # This is normal distribution
#' Anderson_Darling_test(x)
#' #> $is_normality
#' #> [1] TRUE
#' #> $statistic
#' #> [1] 0.1602067
#' #> $pval
#' #> [1] 0.9470751
#'
#' #' x <- stats::runif(100)  # This is non-normal distribution
#' Anderson_Darling_test(x)
#' #> $is_normality
#' #> [1] FALSE
#' #> $statistic
#' #> [1] 1.14719
#' #> $pval
#' #> [1] 0.0050856
#' @references
#' D’Agostino, RalphB., 2017.
#' Goodness-of-Fit Techniques.
#' 1st ed. Routledge. pg. 127 (Table 4.9).
#' https://doi.org/10.1201/9780203753064
Anderson_Darling_test <- function(x){
    x <- sort(x[stats::complete.cases(x)], decreasing = FALSE)
    N <- length(x)

    if (N < 8) stop("sample size must be greater than 7")

    z_score <- (x - mean(x)) / stats::sd(x)

    left_tail <- stats::pnorm(z_score, log.p = TRUE)
    right_tail <- log(1 - stats::pnorm(rev(z_score)))

    # AD: Anderson-Darling statistic
    AD <- -N - mean( (2 * seq_along(x) - 1) * (left_tail + right_tail) )
    statistic <- AD
    AD <- AD * ( 1 + ( 0.75 / N ) + (2.25 / (N ^ 2)) )

    # Piece-wise function
    if (AD <= 0.2)
        pval <- 1 - exp( -13.436 + (101.14 * AD) - (223.73 * AD * AD) )

    if (AD > 0.2 & AD < 0.34)
        pval <- 1 - exp( -8.138 + (42.796 * AD) - (59.938 * AD * AD) )

    if (AD >= 0.34 & AD < 0.6)
        pval <- exp( 0.9177 - (4.279 * AD) - (1.38 * AD * AD) )

    if (AD >= 0.6)
        pval <- exp( 1.2937 - (5.709 * AD) + (0.0186 * AD * AD) )

    if (AD >= 10)
        pval <- 9.9e-12


    is_normality <- pval > 0.05

    ret <- list(
        "is_normality" = is_normality,
        "statistic" = statistic,
        "pval" = pval
    )

    return(ret)
}


Cramer_von_Mises_test <- function(x){
    cat("")
}


D.Angostino_Pearson_test <- function(x){
    cat("")
}


Shapiro_Francia_test <- function(x){
    cat("")
}


Shapiro_Wilk_test <- function(x, method = "original"){
    method <- tolower(method)
    method <- match.arg(method, c("original", "royston"))
}


Jarque_Bera_test <- function(x){
    cat("")
}


Kolmogoro_Smirnov_test <- function(x){
    cat("")
}


Lilliefors_test <- function(x){
    cat("")
}


Pearson_chisquare_test <- function(x){
    cat("")
}


Ryan_Joiner_test <- function(x){
    cat("")
}











