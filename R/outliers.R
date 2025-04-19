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
#' If this is specified, then one-tailed test is performed.
#' If `NULL`, then performed two-tailed test.
#' @param alpha Default: 0.05
#' @param iteration An integer indicating the maximum number of outliers to be detected.
#' Negative values search for all possible outliers (default).
#' @param prob The maximum proportion of the outliers in the dataset. (Default: 0.8)
#'
#' @return A boolean named vector showing the possible outliers.
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- c(round(rnorm(5), 3), -6, 5, 5, round(rnorm(5), 3))
#' Grubbs_test(x)
#' #> -0.626  0.184 -0.836  1.595   0.33     -6      5      5  -0.82  0.487  0.738  0.576 -0.305
#' #>  FALSE  FALSE  FALSE   TRUE  FALSE   TRUE   TRUE   TRUE  FALSE  FALSE  FALSE  FALSE  FALSE
#' #> attr(,"iterations")
#' #> [1] 3
#' #> attr(,"suspect")
#' #> [1] -6.000  5.000  5.000  1.595
#' #> attr(,"G")
#' #> [1] 3.999678 6.237262 2.609002
#' #> attr(,"G_crit")
#' #> [1] 2.462033 2.411560 2.289954
Grubbs_test <- function(
        x,
        xs = NULL,
        alpha = 0.05,
        iteration = -1L,
        prob = 0.8
){
    iteration <- ceiling(iteration)
    if (prob < 0 | prob > 1) stop("The value of `prob` should be within 0 - 1")

    if ( ! is.null(xs) )
        return(.GrubbsTest(x, xs, alpha))

    is_outlier <- logical(length(x))  # A vector contains `FALSE`
    who <- c()
    G <- c()
    G_crit <- c()

    x0 <- x
    i <- 0

    if (iteration <= 0)
        repeat {
            x0 <- x0[ ! is_outlier ]

            out <- .GrubbsTest(x0, xs, alpha)

            is_outlier <- c(unname(out))
            if ( ! any(is_outlier) ) break
            if ( length(x0) < (length(x) * prob) ) break

            i <- i + 1
            who <- append(who, x0[out])
            G <- append(G, attr(out, "G"))
            G_crit <- append(G_crit, attr(out, "G_crit"))
        }

    if (iteration > 0)
        while (i < iteration) {
            x0 <- x0[ ! is_outlier ]

            out <- .GrubbsTest(x0, xs, alpha)

            is_outlier <- c(unname(out))
            if ( ! any(is_outlier) ) break

            i <- i + 1
            who <- append(who, x0[out])
            G <- append(G, attr(out, "G"))
            G_crit <- append(G_crit, attr(out, "G_crit"))
        }

    ret <- structure(
        .Data = ( x %in% who ),
        names = as.character(x),
        iterations = i,
        suspect = who,
        G = G,
        G_crit = G_crit
    )

    return(ret)

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Testing
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    set.seed(1)
    x <- c(round(rnorm(5), 3), -6, 5, 5, round(rnorm(5), 3))
    Grubbs_test(x)
    Grubbs_test(x, xs = -6)
    Grubbs_test(x, iteration = 2)
}


.GrubbsTest <- function(x, xs = NULL, alpha = 0.05)
{
    if ( ! all(is_real_number(x)) ) stop("All input should be real numbers.")
    if ( length(x) < 8 ) warning("Valid observations should be more than 7.")

    x0 <- x

    is_outlier <- logical(length(x0))

    N <- length(x0)
    avg <- mean(x0)
    dif <- abs(x0 - avg)

    two_tail <- is.null(xs)

    alpha <- ifelse(two_tail, alpha / (N + N), alpha / N)

    if (two_tail) {
        suspect_i <- which(dif == max(dif))
        suspect <- x0[suspect_i]

        trim_x <- setdiff(x0, suspect)
        trim_avg <- mean(trim_x)
        trim_std <- stats::sd(trim_x)
        G <- abs(suspect[1] - trim_avg) / trim_std
    } else {
        ## One-tail
        if ( ! is_real_number(xs) || length(xs) > 1 )
            stop("`xs` accepts only one number.")

        if ( ! any(x0 == xs) )
            stop(paste0(xs, " was not found in `x`"))

        suspect_i <- which(dif >= abs(xs - avg))
        suspect <- x0[suspect_i]

        trim_x <- setdiff(x0, suspect)
        trim_avg <- mean(trim_x)
        trim_std <- stats::sd(trim_x)
        G <- abs(xs - trim_avg) / trim_std
    }

    t_crit <- stats::qt(
        p = alpha,
        df = (N - 2),
        lower.tail = FALSE
    )

    G_crit <- ( (N - 1) * t_crit ) / sqrt( N * (N - 2 + t_crit ^ 2) )

    is_outlier[suspect_i] <- ( G > G_crit )

    ret <- structure(
        .Data = is_outlier,
        names = as.character(x0),
        suspect = suspect,
        index = suspect_i,
        G = G,
        G_crit = G_crit
    )

    return(ret)

    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # Testing
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    set.seed(1)
    x <- c(round(rnorm(7, 0, 1), 3), -5, -4, 3, 5, round(rnorm(3, 0, 1), 3))
    .GrubbsTest(x)
    .GrubbsTest(x, xs = 0.738)
}



ROUT_test <- function(x)
{
    cat("Not yet")
}
