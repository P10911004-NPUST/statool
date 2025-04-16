#' Title
#'
#' @param x A numeric vector
#' @param method Detect the outliers with either standard deviation (sd), interquantile range (iqr), or mean (mean)
#' @param n.iqr Default: 1.5
#' @param n.sd Default: 3
#' @param trim Default: 0.1
#' @param lower_quantile Default: 0.2
#' @param upper_quantile Default: 0.8
#' @param use_median Default: TRUE
#'
#' @return A boolean vector
#' @export
#'
#' @examples
#' set.seed(1)
#' x1 <- c(rnorm(5), 4)
#' is_outlier(x1, method = "sd", n.sd = 2)
#' #> FALSE FALSE FALSE FALSE FALSE  TRUE
is_outlier <- function(
        x,
        method = "iqr",
        n.iqr = 1.5,
        n.sd = 3,
        trim = 0.1,
        lower_quantile = 0.2,
        upper_quantile = 0.8,
        use_median = TRUE
){

    if (!is.null(dim(x))) stop("\nInput should be a vector\n")

    method <- match.arg(method, c("sd", "iqr", "mean"))

    x <- stats::na.omit(x)
    stdev <- stats::sd(x)
    center <- mean(x)
    if (use_median) center <- stats::median(x)

    .detect_with_mean <- function(x){
        avg <- mean(x, trim = trim)
        ut <- stats::quantile(x, probs = 1 - (trim / 2), na.rm = TRUE)
        lt <- stats::quantile(x, probs = trim / 2, na.rm = TRUE)
        res <- ( x < lt ) | ( x > ut )
        return(res)
    }

    .detect_with_sd <- function(x){
        ut <- center + n.sd * stdev
        lt <- center - n.sd * stdev
        res <- ( x < lt ) | ( x > ut )
        return(res)
    }

    .detect_with_iqr <- function(x){
        lt <- stats::quantile(x, lower_quantile)
        ut <- stats::quantile(x, upper_quantile)
        IQR <- ut - lt
        res <- (x < (lt - n.iqr * IQR)) | (x > (ut + n.iqr * IQR))
        return(res)
    }

    if (tolower(method) == "mean") res <- .detect_with_mean(x)
    if (tolower(method) == "sd") res <- .detect_with_sd(x)
    if (tolower(method) == "iqr") res <- .detect_with_iqr(x)

    return(res)
}



#' Grubbs's test
#'
#' Description
#'
#' @import stats
#' @param x A numeric vector (normal distribution) to test for an outlier.
#' @param xs Which observation is suspected as an outlier.
#' If this is specified, then one-tailed test is performed. If NULL, then performed two-tailed test.
#' @param alpha Default: 0.05
#'
#' @return A list containing 3 vectors:
#' 1. G: The statistic value for the test.
#' 2. Gcrit: The threshold of the test to reject null hypothesis.
#' 3. is_outlier: Boolean values indicating which one is the outlier (G > Gcrit).
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- c(round(rnorm(10), 2), 5)
#' Grubbs_test(x)
#' #> $G
#' #> [1] 2.689996
#' #>
#' #> $Gcrit
#' #> [1] 2.35473
#' #>
#' #> $is_outlier
#' #> -0.63  0.18 -0.84   1.6  0.33 -0.82  0.49  0.74  0.58 -0.31     5
#' #> FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE
Grubbs_test <- function(x, xs = NULL, alpha = 0.05)
{
    x0 <- stats::na.omit(x)

    N <- length(x0)
    avg <- mean(x0)
    std <- stats::sd(x0)

    two_tail <- is.null(xs)
    alpha <- ifelse(two_tail, alpha / (2 * N), alpha / N)

    d <- ifelse(two_tail, abs(x0 - avg), abs(xs - avg))

    if (two_tail) {
        d <- abs(x0 - avg)
        max_i <- which.max(d)
        G <- d[max_i] / std
    } else {
        G <- abs(xs - avg) / std
    }

    t_crit <- stats::qt(
        p = alpha,
        df = N - 2,
        lower.tail = FALSE
    )

    G_crit <- ( (N - 1) * t_crit ) / sqrt( N * (N - 2 + t_crit * t_crit) )

    if (two_tail) {
        is_outlier <- stats::setNames(logical(length(x0)), as.character(x0))
        is_outlier[max_i] <- (G > G_crit)
    } else {
        is_outlier <- (G > G_crit)
    }

    x1 <- x0[ ! is_outlier ]
    if ( ! is_normality(x1) )
        warning("Not normally distributed. Grubbs's test may be NOT suitable for this data.")

    ret <- list(
        G = G,
        # tcrit = t_crit,
        Gcrit = G_crit,
        is_outlier = is_outlier
    )

    return(ret)

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Testing ====
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    set.seed(1)
    x <- c(round(rnorm(10), 2), 5)
    Grubbs_test(x)
}



ROUT_test <- function(x)
{
    cat("Not yet")
}
