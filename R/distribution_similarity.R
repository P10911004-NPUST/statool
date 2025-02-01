# Test whether the samples come from the same distribution

# k-sample Anderson Darling test ====
ksample_anderson_darling <- function(data, formula, midrank = TRUE)
{
    df0 <- stats::model.frame(formula = formula, data = data)
    colnames(df0) <- c("y", "x")

    df0 <- df0[stats::complete.cases(df0), ]
    df0 <- df0[order(df0$y, decreasing = FALSE), ]

    # n: number of samples of each groups
    group_names <- unique(df0$x)
    group_sizes <- tapply(df0$y, df0$x, "length")
    k <- length(group_names)  # number of groups
    N <- nrow(df0)  # total number of observations

    Z <- sort(df0$y, decreasing = FALSE)
    Zstar <- sort(unique(Z), decreasing = FALSE)

    if (length(Zstar) < 2) stop("Needs more than one distinct observation")
    if (k < 2) stop("Require at least 2 groups")

    if (midrank) {
        i <- tapply(df0$y, df0$x, function(x) .anderson_midrank(x, Z, Zstar, N))
        i_cumsum <- utils::tail(base::cumsum(i), 1)
        AakN2 <- i_cumsum * ((N - 1) / N)
    } else {
        cat("Not yet")
    }

    AkN2 <- AakN2



    H_ <- base::sum(1 / group_sizes)
    h_ <- base::sum(1 / seq(1, N - 1, 1))
    g_ <- NULL
    a_ <- NULL
    b_ <- NULL
    c_ <- NULL
    d_ <- NULL
    # sigmaN2 <- (a * (N ^ 3)) / ()

    return(AakN2)
}


## Equation 6 - EDF (version 1) ====
.anderson_right_EDF <- function(x, Z, Zstar, N)
{
    cat("Not yet")
}


## Equation 7 - Midrank (version 2) ====
.anderson_midrank <- function(x, Z, Zstar, N)
{
    x <- sort(x, decreasing = FALSE)
    n <- length(x)  # number of observations for each group

    x_left <- search_sorted(Zstar, x, side = "left")
    x_right <- search_sorted(Zstar, x, side = "right")

    Z_left <- search_sorted(Zstar, Z, side = "left")
    Z_right <- search_sorted(Zstar, Z, side = "right")  ## have problem

    lj <- 1
    if (N != length(Zstar)) lj <- Z_right - Z_left  # If tied-values exist

    Baj <- Z_left + (lj / 2)
    fij <- x_right - x_left
    Maij <- x_right - (fij / 2)

    j_sum <- (lj / N) * ( (N * Maij - n * Baj) ^ 2 ) / ( (Baj * (N - Baj)) - (N * lj / 4) )
    j_sum <- sum(j_sum)

    i <- j_sum / n

    return(i)
}




# ## sigma_N for k-sample Anderson Darling ====
# .sigmaN <- function(){
#
# }



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Testing ====
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
source("./numpy.R")
exp_data <- data.frame(
    A = c(38.7, 41.5, 43.8, 44.5, 45.5, 46.0, 47.7, 58.0),
    B = c(39.2, 39.3, 39.7, 41.4, 41.8, 42.9, 43.3, 45.8),
    C = c(34.0, 35.0, 39.0, 40.0, 43.0, 43.0, 44.0, 45.0),
    D = c(34.0, 34.8, 34.8, 35.4, 37.2, 37.8, 41.2, 42.8)
)

res0 <- with(exp_data, kSamples::ad.test(A, B, C, D))
res0
res0$ad[2, ][["AD"]]

exp_data <- stats::reshape(
    data = exp_data,
    direction = "long",
    timevar = "Laboratory",
    times = colnames(exp_data),
    varying = colnames(exp_data),
    v.names = "Smoothness"
)

x <- 0
Z <- c(
    34.0, 34.0, 34.8, 34.8, 35.0, 35.4, 37.2, 37.8, 38.7, 39.0, 39.2, 39.3,
    39.7, 40.0, 41.2, 41.4, 41.5, 41.8, 42.8, 42.9, 43.0, 43.0, 43.3, 43.8,
    44.0, 44.5, 45.0, 45.5, 45.8, 46.0, 47.7, 58.0
)
Zstar <- c(
    34.0, 34.8, 35.0, 35.4, 37.2, 37.8, 38.7, 39.0, 39.2, 39.3, 39.7, 40.0,
    41.2, 41.4, 41.5, 41.8, 42.8, 42.9, 43.0, 43.3, 43.8, 44.0, 44.5, 45.0,
    45.5, 45.8, 46.0, 47.7, 58.0
)
N <- 32

# tapply(exp_data$Smoothness, exp_data$Laboratory, function(x) .anderson_midrank(x, Z, Zstar, N))

res0 <- ksample_anderson_darling(exp_data, Smoothness ~ Laboratory)
res0











