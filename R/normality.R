# Normality test ====


#' Skewness of the sample / population distribution
#'
#' Return the skewness of the sample / population distribution.
#' The skewness value for a normal distribution should be close to zero.
#' A negative skewness value indicates left-skewed distribution (concentrated on the right);
#' A positive skewness value indicates right-skewed distribution (concentrated on the left).
#'
#' @param x A vector of numeric values. NA values will be discarded automatically.
#' @param .population Calculate the skewness of the population or the sample.
#'
#' @return Numeric value
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- stats::rnorm(10, 10)
#' skewness(x)
#' #> 0.3512426
skewness <- function(x, .population = FALSE){
    x <- x[stats::complete.cases(x)]
    N <- length(x)
    if (N < 3 || is.null(N)) stop("Sample size should be greater than 3")

    if (.population){
        ret <- sum( ( (x - mean(x)) / sd_population(x) ) ^ 3 ) / N  # sd_population() <- utils.R
    } else {
        numerator <- N * sum( (x - mean(x)) ^ 3 )
        denominator <- (N - 1) * (N - 2) * ( stats::sd(x) ^ 3 )
        ret <- numerator / denominator
    }

    return(ret)
}



#' Kurtosis of the sample / population distribution
#'
#' Return the kurtosis of the sample / population distribution.
#'
#' @param x A vector of numeric values. NA values will be discarded automatically.
#' @param .population Calculate the kurtosis of the population or the sample.
#'
#' @return Numeric value
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- stats::rnorm(10, 10)
#' kurtosis(x)
#' #> -0.3169031
kurtosis <- function(x, .population = FALSE){
    x <- x[stats::complete.cases(x)]
    N <- length(x)

    if (.population) {
        ret <- sum( ( (x - mean(x)) / sd_population(x) ) ^ 4 ) / N
    } else {
        # Block A >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        numerator <- (N + 1) * N * sum( (x - mean(x)) ^ 4 )
        denominator <- (N - 1) * (N - 2) * (N - 3) * ( stats::sd(x) ^ 4 )
        block_A <- numerator / denominator

        # Block B >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        numerator <- 3 * ( (N - 1) ^ 2 )
        denominator <- (N - 2) * (N - 3)
        block_B <- numerator / denominator

        ret <- block_A - block_B
    }

    return(ret)
}



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Empirical Distribution Function (EDF) tests ====
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

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
    cat("Not yet")
}


Kolmogoro_Smirnov_test <- function(x){
    cat("Not yet")
}


Lilliefors_test <- function(x){
    cat("Not yet")
}




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Regression and Correlation tests ====
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Shapiro_Wilk_test <- function(x){
    cat("Not yet")
}


Shapiro_Wilk_Royston_test <- function(x){
    cat("Not yet")
}


Shapiro_Francia_test <- function(x){
    cat("Not yet")
}


Ryan_Joiner_test <- function(x){
    cat("Not yet")
}




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Moment tests ====
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
D.Angostino_Pearson_test <- function(x){
    cat("Not yet")
}


Jarque_Bera_test <- function(x){
    cat("Not yet")
}




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Other tests ====
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Pearson_chisquare_test <- function(x){
    cat("Not yet")
}
