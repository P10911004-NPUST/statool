# Feature scaling (data normalization)

min_max_scaling <- function(x, range = c(0, 1)){
    if (!is.vector(x))
        stop("Input `x` should be a numeric vector.\n")
    if (length(x) < 3)
        stop("Length of `x` should be > 2.\n")
    if (!is.vector(range) || length(range) != 2)
        stop("`range` should be a numeric vector with length of two.\n")

    x <- x[stats::complete.cases(x)]
    x_min <- min(x)
    x_max <- max(x)

    ret <- ((x - x_min) / (x_max - x_min)) * (range[2] - range[1]) + range[1]

    return(ret)
}


mean_scaling <- function(x, range = c(0, 1)){
    if (!is.vector(x))
        stop("Input `x` should be a numeric vector.\n")
    if (length(x) < 3)
        stop("Length of `x` should be > 2.\n")
    if (!is.vector(range) || length(range) != 2)
        stop("`range` should be a numeric vector with length of two.\n")

    x <- x[stats::complete.cases(x)]
    x_min <- min(x)
    x_max <- max(x)

    ret <- ((x - mean(x)) / (x_max - x_min)) * (range[2] - range[1]) + range[1]

    return(ret)
}


zscore_scaling <- function(x){
    if (!is.vector(x))
        stop("Input `x` should be a numeric vector.\n")
    if (length(x) < 3)
        stop("Length of `x` should be > 2.\n")

    x <- x[stats::complete.cases(x)]
    z <- (x - mean(x)) / stats::sd(x)

    return(z)
}


robust_scaling <- function(x){
    if (!is.vector(x))
        stop("Input `x` should be a numeric vector.\n")
    if (length(x) < 5)
        stop("Length of `x` should be > 4.\n")

    x <- x[stats::complete.cases(x)]

    ret <- (x - stats::median(x)) / stats::IQR(x)

    return(ret)
}


unit_vector_scaling <- function(x, method = "Lp"){
    method <- match.arg(method, c("L1", "L2", "Lp"))
    if (!is.vector(x))
        stop("Input `x` should be a numeric vector.\n")
    if (length(x) < 5)
        stop("Length of `x` should be > 4.\n")

    x <- x[stats::complete.cases(x)]

    cat("Not yet")
}








