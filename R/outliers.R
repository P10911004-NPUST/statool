is_outlier <- function(
        x,
        method = "iqr",
        n.iqr = 1.5,
        n.sd = 3,
        trim = 0.1,
        lower_quantile = 0.2,
        upper_quantile = 0.8,
        use_median = FALSE
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
