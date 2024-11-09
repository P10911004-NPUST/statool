# Feature scaling (data normalization)

#' Minimum-Maximum (Min-Max) scaling
#'
#' Feature scaling method normalizes the input data to a desired range (generally from 0 to 1).
#' The transformation was done by:
#' `[(x - min(x)) / (max(x) - min(x))] * (range_ceiling - range_floor) + range_floor`
#'
#' @param x A numeric vector.
#' @param range A numeric vector with length of 2. The desired range for the data normalization.
#' The default is c(0, 1).
#'
#' @return A numeric vector, i.e., the normalized data.
#' @export
#' @author Joon-Keat Lai
#'
#' @examples
#' set.seed(1)
#' x <- stats::rnorm(10)
#' min_max_scaling(x)
#' #> [1] 0.086047962 0.419296552 0.000000000 1.000000000 0.479300618
#' #> [6] 0.006236443 0.544264487 0.647475101 0.580609856 0.218124222
#'
#' min_max_scaling(x, range = c(-2, 2))
#' #> [1] -1.65580815 -0.32281379 -2.00000000  2.00000000 -0.08279753
#' #> [6] -1.97505423  0.17705795  0.58990040  0.32243942 -1.12750311
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


#' Mean scaling
#'
#' Feature scaling method normalizes the input data to a desired range (generally from 0 to 1).
#' The transformation was done by:
#' `[(x - mean(x)) / (max(x) - min(x))] * (range_ceiling - range_floor) + range_floor`
#'
#' @param x A numeric vector.
#' @param range A numeric vector with length of 2. The desired range for the data normalization.
#' The default is c(0, 1).
#'
#' @return A numeric vector, i.e., the normalized data.
#' @export
#' @author Joon-Keat Lai
#'
#' @examples
#' set.seed(1)
#' x <- stats::rnorm(10)
#' mean_scaling(x)
#' #> [1] -0.31208756  0.02116103 -0.39813552  0.60186448  0.08116509
#' #> [6] -0.39189908  0.14612896  0.24933958  0.18247433 -0.18001130
#'
#' mean_scaling(x, range = c(-2, 2))
#' #> [1] -3.2483502 -1.9153559 -3.5925421  0.4074579 -1.6753396
#' #> [6] -3.5675963 -1.4154841 -1.0026417 -1.2701027 -2.7200452
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


#' Z-score scaling (Standardization)
#'
#' Feature scaling method standardizes the input data to Z-score.
#' The transformation was done by:
#' `(x - mean(x)) / sd(x)`
#'
#' @param x A numeric vector.
#'
#' @return A numeric vector, i.e., the standardized data.
#' @export
#' @author Joon-Keat Lai
#'
#' @examples
#' set.seed(1)
#' x <- stats::rnorm(10)
#' z_score_scaling(x)
#' #> [1] -0.97190653  0.06589991 -1.23987805  1.87433300  0.25276523
#' #> [6] -1.22045645  0.45507643  0.77649606  0.56826358 -0.56059319
z_score_scaling <- function(x){
    if (!is.vector(x))
        stop("Input `x` should be a numeric vector.\n")
    if (length(x) < 3)
        stop("Length of `x` should be > 2.\n")

    x <- x[stats::complete.cases(x)]
    z <- (x - mean(x)) / stats::sd(x)

    return(z)
}



#' Robust scaling
#'
#' Feature scaling method normalizes the input data by the
#' median and the inter-quantile range (IQR, P75 - P25).
#' The transformation was done by:
#' `(x - median(x)) / IQR(x)`
#'
#' @param x A numeric vector.
#'
#' @return A numeric vector, i.e., the normalized data.
#' @export
#' @author Joon-Keat Lai
#'
#' @examples
#' set.seed(1)
#' x <- stats::rnorm(10)
#' robust_scaling(x)
#' #> [1] -0.80284101 -0.06630921 -0.99302054  1.21713675  0.06630921
#' #> [6] -0.97923702  0.20988958  0.43800127  0.29021856 -0.51093170
robust_scaling <- function(x){
    if (!is.vector(x))
        stop("Input `x` should be a numeric vector.\n")
    if (length(x) < 5)
        stop("Length of `x` should be at least 5.\n")

    x <- x[stats::complete.cases(x)]

    ret <- (x - stats::median(x)) / stats::IQR(x)

    return(ret)
}


unit_vector_scaling <- function(x, method = "Lp"){
    method <- toupper(method)
    method <- match.arg(method, c("L1", "L2", "LP"))
    if (!is.vector(x))
        stop("Input `x` should be a numeric vector.\n")
    if (length(x) < 5)
        stop("Length of `x` should be > 4.\n")

    x <- x[stats::complete.cases(x)]

    if (method == "L1")
        ret <- x / sum(abs(x))

    if (method == "L2")
        ret <- x / sqrt(sum(x ^ 2))

    if (method == "LP")
        cat("Not yet")

    return(ret)
}








