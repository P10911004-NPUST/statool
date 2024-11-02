# Normality test ====


#' Skewness of the sample / population distribution
#'
#' Return the skewness of the sample / population distribution.
#' The skewness of a normal distribution approximately equal to zero.
#' A negative skewness value indicates left-skewed distribution (concentrated on the right);
#' A positive skewness value indicates right-skewed distribution (concentrated on the left).
#'
#' @param x A vector of numeric values. NA values will be discarded automatically.
#' @param .population Calculate the skewness of the population or the sample.
#'
#' @return Numeric value
#'
#' @export
#' @author Joon-Keat Lai
#'
#' @examples
#' set.seed(1)
#' x <- stats::rnorm(10, 10)
#' skewness(x)
#' #> 0.3512426
skewness <- function(x, .population = FALSE) {
    if ( !is_vector(x) || !is.numeric(x) ) stop("x should be a numeric vector")

    x <- x[stats::complete.cases(x)]
    N <- length(x)
    if (N < 3 || is.null(N)) stop("Sample size must be greater than 3")

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
#' Normally distributed data has kurtosis around 3.0.
#'
#' @param x A vector of numeric values. NA values will be discarded automatically.
#' @param excess_kurtosis Logical (default: TRUE). Calculate the excess kurtosis (Ke); Ke = K - 3.
#' @param .population Logical (default: FALSE). Calculate the kurtosis of the population or the sample.
#'
#' @return Numeric value
#'
#' @export
#' @author Joon-Keat Lai
#'
#' @examples
#' set.seed(1)
#' x <- stats::rnorm(10, 10)
#' kurtosis(x)
#' #> -0.3169031
kurtosis <- function(x, excess_kurtosis = TRUE, .population = FALSE) {
    if ( !is_vector(x) || !is.numeric(x) ) stop("x should be a numeric vector")

    x <- x[stats::complete.cases(x)]
    N <- length(x)
    if (N < 3 || is.null(N)) stop("Sample size must be greater than 3")

    if (.population) {
        ret <- sum( ( (x - mean(x)) / sd_population(x) ) ^ 4 ) / N
    } else {
        # Block A >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        A_numerator <- (N + 1) * N * sum( (x - mean(x)) ^ 4 )
        A_denominator <- (N - 1) * (N - 2) * (N - 3) * ( stats::sd(x) ^ 4 )
        block_A <- A_numerator / A_denominator

        # Block B >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        B_numerator <- 3 * ( (N - 1) ^ 2 )
        B_denominator <- (N - 2) * (N - 3)
        block_B <- B_numerator / B_denominator

        ret <- block_A - block_B
    }

    if ( isFALSE(excess_kurtosis) ) ret <- ret + 3

    return(ret)
}



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Empirical Distribution Function (EDF) tests ====
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


#' Anderson-Darling normality test
#'
#' The Anderson-Darling test is an empirical distribution function (EDF) omnibus test
#' for the hypothesis of normality. The test statistic is
#' \deqn{
#' A = -n -\frac{1}{n} \sum_{i=1}^{n} [2i-1]
#' [\ln(p_{(i)}) + \ln(1 - p_{(n-i+1)})],
#' }
#' where
#' \eqn{p_{(i)} = \Phi([x_{(i)} - \overline{x}]/s)}.
#' \eqn{\Phi} is the cumulative distribution function of the standard normal distribution;
#' \eqn{\overline{x}} and \eqn{s} are mean and standard deviation of the data values.
#' The p-value is computed from the modified statistic
#' \eqn{Z=A (1.0 + 0.75/n +2.25/n^{2})} according to the Table 4.9 in D’Agostino (2017).
#'
#' @param x A numeric vector
#'
#' @return A list with three elements:
#' 1. is_normality: logical value. TRUE indicates x is a normal distribution.
#' 2. statistic: the Anderson-Darling test statistic value.
#' The higher the value, the lower the probability of normality.
#' 3. pval: the p-value for the test. pval > 0.05 indicates normal distribution.
#'
#' @export
#' @author Joon-Keat Lai
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
#'
#' @references
#' D’Agostino, R.B. 2017.
#' Goodness-of-Fit Techniques.
#' 1st ed. Routledge. pg. 127 (Table 4.9).
#' https://doi.org/10.1201/9780203753064
Anderson_Darling_test <- function(x) {
    if ( !is_vector(x) || !is.numeric(x) ) stop("x should be a numeric vector")

    x <- sort(x[stats::complete.cases(x)], decreasing = FALSE)
    N <- length(x)

    if ( N < 8 || is.null(N) ) stop("Sample size must be greater than 7")

    z_score <- (x - mean(x)) / stats::sd(x)

    left_tail <- stats::pnorm(z_score, log.p = TRUE)
    right_tail <- log(1 - stats::pnorm(rev(z_score)))

    # AD: Anderson-Darling test statistic
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

    is_normality <- pval >= 0.05

    ret <- list(
        "is_normality" = is_normality,
        "statistic" = statistic,
        "pval" = pval
    )

    return(ret)
}



#' Cramer-von Mises normality test
#'
#' The Cramer-von Mises test is an empirical distribution function (EDF) omnibus test
#' for the hypothesis of normality. The test statistic is
#' \deqn{
#' W = \frac{1}{12 n} + \sum_{i=1}^{n} \left(p_{(i)} - \frac{2i-1}{2n}\right)^2,
#' }{W = 1/(12n)  +  \sum_{i=1}^n (p_(i) - (2i-1)/(2n))^2,} \
#' where \eqn{p_{(i)} = \Phi([x_{(i)} - \overline{x}]/s)},
#' \eqn{\Phi} is the cumulative distribution function of the standard normal distribution,
#' \eqn{\overline{x}} and \eqn{s} are the mean and standard deviation.
#' The p-value is computed from the modified statistic \eqn{Z = W (1.0 + 0.5/n)}
#' according to Table 4.9 in Stephens (1986).
#'
#' @param x A numeric vector
#'
#' @return A list with three elements:
#' 1. is_normality: logical value. TRUE indicates x is a normal distribution.
#' 2. statistic: the Cramer-von Mises test statistic value (W).
#' The higher the value, the lower the probability of normality.
#' 3. pval: the p-value for the test. pval > 0.05 indicates normal distribution.
#'
#' @export
#' @author Joon-Keat Lai
#'
#' @examples
#' set.seed(1)
#' x <- stats::rnorm(100)  # This is normal distribution
#' Cramer_von_Mises_test(x)
#' #> $is_normality
#' #> [1] TRUE
#' #> $statistic
#' #>         W
#' #> 0.0260313
#' #> $pval
#' #> [1] 0.8945191
#'
#' x <- stats::runif(100)  # This is non-normal distribution
#' Cramer_von_Mises_test(x)
#' #> $is_normality
#' #> [1] FALSE
#' #> $statistic
#' #>         W
#' #> 0.2746561
#' #> $pval
#' #> [1] 0.0006342302
#'
#' @references
#' Stephens, M.A. 1986.
#' Tests based on EDF statistics.
#' Eds.: D'Agostino, R.B. and Stephens, M.A.
#' Goodness-of-Fit Techniques.
#' Marcel Dekker, New York.
#'
#' Thode Jr., H.C. 2002.
#' Testing for  Normality.
#' Marcel Dekker, New York.
Cramer_von_Mises_test <- function(x) {
    x <- sort(x[stats::complete.cases(x)])
    N <- length(x)
    i <- seq_along(x)

    if (N <= 7) stop("sample size must be greater than 7")

    z <- (x - mean(x)) / stats::sd(x)
    p <- stats::pnorm(z)

    block_A <- 1 / (12 * N)
    block_B <- (2 * i - 1) / (2 * N)

    W <- block_A + sum( (p - block_B) ^ 2 )
    statistic <- c("W" = W)

    W <- (1 + 0.5 / N) * W

    if (W >= 1.1)
        pval <- 9.9e-10

    if (W < 1.1 && W >= 0.092)
        pval <- exp(1.111 - 34.242 * W + 12.832 * W ^ 2)

    if (W < 0.092 && W >= 0.051)
        pval <- exp(0.886 - 31.62 * W + 10.897 * W ^ 2)

    if (W < 0.051 && W >= 0.0275)
        pval <- 1 - exp(-5.903 + 179.546 * W - 1515.29 * W ^ 2)

    if (W < 0.0275)
        pval <- 1 - exp(-13.953 + 775.5 * W - 12542.61 * W ^ 2)

    is_normality <- pval >= 0.05

    ret <- list(
        "is_normality" = is_normality,
        "statistic" = statistic,
        "pval" = pval
    )

    return(ret)
}


Kolmogoro_Smirnov_test <- function(x) {
    cat("Not yet")
}


Lilliefors_test <- function(x) {
    cat("Not yet")
}




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Regression and Correlation tests ====
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Shapiro_Wilk_test <- function(x) {
    cat("Not yet")
}


Shapiro_Wilk_Royston_test <- function(x) {
    cat("Not yet")
}



#' Shapiro-Francia normality test
#'
#' The Shapiro-Francia test is a type of regression and correlation test for normality
#' evaluation. It estimates the squared correlation between the ordered sample values
#' and the approximate expected ordered quantiles from the standard normal distribution.
#' The p-value is computed from the formula given by Royston (1993).
#'
#' @param x A numeric vector
#'
#' @return A list with three elements:
#' 1. is_normality: logical value. TRUE indicates x is a normal distribution.
#' 2. statistic: the Shapiro-Francia test statistic value (W).
#' The higher the value, the lower the probability of normality.
#' 3. pval: the p-value for the test. pval > 0.05 indicates normal distribution.
#'
#' @export
#' @author Joon-Keat Lai
#'
#' @examples
#' set.seed(1)
#' x <- stats::rnorm(100)  # This is normal distribution
#' Shapiro_Francia_test(x)
#' #> $is_normality
#' #> [1] TRUE
#' #> $statistic
#' #>         W
#' #> 0.9958174
#' #> $pval
#' #> [1] 0.9718829
#'
#' x <- stats::runif(100)  # This is non-normal distribution
#' Shapiro_Francia_test(x)
#' #> $is_normality
#' #> [1] FALSE
#' #> $statistic
#' #>         W
#' #> 0.9386237
#' #> $pval
#' #> [1] 0.0003371829
#'
#' @references
#' Royston, P. 1993.
#' A pocket-calculator algorithm for the Shapiro-Francia test for non-normality: an application to medicine.
#' Statistics in Medicine. 12, 181-184.
#'
#' Thode Jr., H.C. 2002.
#' Testing for  Normality.
#' Marcel Dekker, New York.
Shapiro_Francia_test <- function(x) {
    x <- sort(x[stats::complete.cases(x)], decreasing = FALSE)
    N <- length(x)

    if ( N < 5 || N > 5000 ) stop("sample size must be between 5 and 5000")

    y <- stats::qnorm( stats::ppoints( N, a = (3 / 8) ) )
    W <- stats::cor(x, y) ^ 2
    u <- log(N)
    v <- log(u)
    mu <- (1.0521 * (v - u)) - 1.2725
    sig <- 1.0308 - 0.26758 * (v + 2 / u)

    z <- ( log(1 - W) - mu ) / sig

    pval <- stats::pnorm(z, lower.tail = FALSE)

    is_normality <- pval >= 0.05

    ret <- list(
        "is_normality" = is_normality,
        "statistic" = c("W" = W),
        "pval" = pval
    )

    return(ret)
}



#' Ryan-Joiner normality test
#'
#' Test the data normality (normal distribution).
#'
#' @param x A numeric vector
#' @param alpha The error tolerance, either 0.1, 0.05 (default), or 0.01.
#'
#' @return A list with three elements
#' 1. is_normality: logical value. TRUE: x is normal distribution
#' 2. statistic: the Ryan-Joiner test statistic (Rp) and its critical value (Rc).
#' Rp > Rc accepts normality assumption.
#' 3. pval: return an `NA`. the p-value of the statistics is currently not available.
#'
#' @export
#' @author Joon-Keat Lai
#'
#' @examples
#' set.seed(1)
#' x <- stats::rnorm(100)  # This is normal distribution
#' Ryan_Joiner_test(x)
#' #> $is_normality
#' #> [1] TRUE
#' #> $statistic
#' #>        Rp           Rc
#' #> 0.9979065    0.9874371
#' #> $pval
#' #> [1] NA
#'
#' set.seed(1)
#' x <- stats::runif(100, 1)  # This is non-normal distribution
#' Ryan_Joiner_test(x)
#' #> $is_normality
#' #> [1] FALSE
#' #> $statistic
#' #>        Rp            Rc
#' #>      -Inf     0.9874371
#' #> $pval
#' #> [1] NA
#'
#' @references
#' Ryan T.A. and Joiner B.L. 1976.
#' Normal probability plots and tests for normality.
#' Technical Report.
#' Statistics Department, The Pennsylvania State University.
Ryan_Joiner_test <- function(x, alpha = 0.05) {
    if ( !alpha %in% c(0.1, 0.05, 0.01) ) stop("`alpha` should be either 0.1, 0.05, 0.01")
    if ( !is_vector(x) || !is.numeric(x) ) stop("`x` should be a numeric vector.")

    x <- x[stats::complete.cases(x)]
    x <- sort(x, decreasing = FALSE)
    N <- length(x)
    zero_index <- seq_along(x) - 1
    ti <- (zero_index + 1 - 0.375) / (N + 0.25)
    bi <- stats::qnorm(ti)

    numerator <- sum(x * bi)
    denominator <- sqrt( sum( (x - mean(x)) ^ 2 ) * sum(bi ^ 2) )
    Rp <- numerator / denominator

    # Critical value for Rp
    if (alpha == 0.10) Rc <- c(b0 = 1.0071, b1 = -0.1271, b2 = -0.3682, b3 = 0.7780)
    if (alpha == 0.05) Rc <- c(b0 = 1.0063, b1 = -0.1288, b2 = -0.6118, b3 = 1.3505)
    if (alpha == 0.01) Rc <- c(b0 = 0.9963, b1 = -0.0211, b2 = -1.4106, b3 = 3.1791)

    Rc <- unname(Rc)
    Rc <- Rc[1] + (Rc[2] / sqrt(N)) + (Rc[3] / N) + (Rc[4] / (N ^ 2))

    statistic <- c("Rp" = Rp, "Rc" = Rc)
    is_normality <- Rp > Rc

    ret <- list(
        "is_normality" = is_normality,
        "statistic" = statistic,
        "pval" = NA
    )

    return(ret)
}


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Moment tests ====
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
D.Agostino_Pearson_test <- function(x) {
    cat("Not yet")
}



#' Jarque-Bera normality test
#'
#' Test the data normality (normal distribution).
#' Insensitive to uniform distribution. Not recommended.
#'
#' @param x A numeric vector
#'
#' @return A list with three elements
#' 1. is_normality: logical value. TRUE: x is normal distribution
#' 2. statistic: the Jarque-Bera test statistic value.
#' The higher the value, the lower the probability of normality.
#' 3. pval: the p-value of the statistics. pval >= 0.05 indicates normal distribution.
#'
#' @export
#' @author Joon-Keat Lai
#'
#' @examples
#' set.seed(1)
#' x <- stats::rnorm(100)  # This is normal distribution
#' Jarque_Bera_test(x)
#' #> $is_normality
#' #> [1] TRUE
#' #> $statistic
#' #> [1] 0.1103686
#' #> $pval
#' #> [1] 0.9463107
#'
#' set.seed(1)
#' x <- stats::rchisq(100, 1)  # This is non-normal distribution
#' Jarque_Bera_test(x)
#' #> $is_normality
#' #> [1] FALSE
#' #> $statistic
#' #> [1] 615.3756
#' #> $pval
#' #> [1] 2.359915e-134
#'
#' @references
#' Jarque, C.M. and Bera, A.K. 1980.
#' Efficient tests for normality, homoscedasticity and serial independence of regression residuals.
#' Economics letters, 6(3), 255-259.
Jarque_Bera_test <- function(x) {
    if ( !is_vector(x) || !is.numeric(x) ) stop("x should be a numeric vector")

    x <- x[stats::complete.cases(x)]
    N <- length(x)
    if ( N < 3 || is.null(N) ) stop("Sample size must be greater than 3")

    skew <- skewness(x)
    Ke <- kurtosis(x, excess_kurtosis = TRUE)

    JB <- (N / 6) * ( (skew ^ 2) + ((Ke ^ 2) / 4) )

    pval <- exp(-JB / 2)  # equivalent to `pchisq(JB, 2, lower.tail = FALSE)`

    is_normality <- pval >= 0.05

    ret <- list(
        "is_normality" = is_normality,
        "statistic" = JB,
        "pval" = pval
    )

    return(ret)
}




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Other tests ====
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Pearson_chisquare_test <- function(x) {
    cat("Not yet")
}
