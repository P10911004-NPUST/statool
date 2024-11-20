# Test whether the samples come from the same distribution

ksample_anderson_darling <- function(data, formula){

    df0 <- stats::model.frame(formula = formula, data = data)
    colnames(df0) <- c("y", "x")

    df0 <- df0[stats::complete.cases(df0), ]
    df0 <- df0[order(df0$y, decreasing = FALSE), ]

    # k: number of groups
    # N: total number of observations
    # n: number of samples of each groups
    group_names <- unique(df0$x)
    k <- length(group_names)
    N <- nrow(df0)

    Z <- sort(df0$y, decreasing = FALSE)
    Zstar <- sort(unique(Z), decreasing = FALSE)

    if (length(Zstar) < 2) stop("Needs more than one distinct observation")
    if (k < 2) stop("Require at least 2 groups")

    stats_mat <- tapply(df0$y, df0$x, function(x) anderson_midrank(x, Z, Zstar, N))
    # Z_left <- search_sorted(Z, Zstar, side = "left")
    # Z_right <- search_sorted(Z, Zstar, side = "right")
    #
    # if (N == length(Zstar)) lj <- 1  # If all observations are unique values
    # if (N != length(Zstar)) lj <- Z_right - Z_left  # If tied-values exist
    #
    # Baj <- Z_left + (lj / 2)

    # inner <- c()
    # for (g in group_names){
    #     ni <- group_sizes[[g]]
    #     xi <- subset(df0, x == g, y)[, , drop = TRUE]
    #     x_right <- search_sorted(xi, Zstar, side = "right")
    #     x_left <- search_sorted(xi, Zstar, side = "left")
    #     fij <- x_right - x_left
    #     Maij <- x_right - (fij / 2)
    #
    #     block_A <- 1 / ni
    #     block_B <- lj / N
    #     block_C <- ( (N * Maij - ni * Baj) ^ 2 ) / ( (Baj * N) - (Baj ^ 2) - (N * lj / 4) )
    #     inner <- append( inner, block_A * sum(block_B * block_C) )
    # }

    # AakN2 <- ((N - 1) / N) * sum(inner)

    # .AakN2 <- function(group_name){
    #     ni <- group_sizes[[group_name]]
    #     xi <- subset(df0, x == group_name, y)[, , drop = TRUE]
    #     x_right <- search_sorted(xi, Zstar, side = "right")
    #     x_left <- search_sorted(xi, Zstar, side = "left")
    #     fij <- x_right - x_left
    #     Maij <- x_right - (fij / 2)
    #
    #     block_A <- 1 / ni
    #     block_B <- lj / N
    #     block_C <- ( (N * Maij - ni * Baj) ^ 2 ) / ( (Baj * N) - (Baj ^ 2) - (N * lj / 4) )
    #     inner <- block_A * sum(block_B * block_C)
    #     return(inner)
    # }
    #
    # ret <- sapply(group_names, .AakN2)


    return(stats_mat)
}


anderson_midrank <- function(x, Z, Zstar, N){
    x <- sort(x, decreasing = FALSE)
    n <- length(x)

    x_left <- search_sorted(x, Zstar, side = "left")
    x_right <- search_sorted(x, Zstar, side = "right")
    Z_left <- search_sorted(Z, Zstar, side = "left")
    Z_right <- search_sorted(Z, Zstar, side = "right")

    lj <- 1
    if (N != length(Zstar)) lj <- Z_right - Z_left  # If tied-values exist

    Baj <- Z_left + (lj / 2)
    fij <- x_right - x_left
    Maij <- x_right - (fij / 2)

    j_sum <- (lj / N) * ( (N * Maij - n * Baj) ^ 2 ) / ( (Baj * (N - Baj)) - (N * lj / 4) )
    j_sum <- sum(j_sum)

    return(j_sum)
}


# u1 <- c(1.0066, -0.9587,  0.3462, -0.2653, -1.3872)
# u2 <- c(0.1005, 0.2252, 0.4810, 0.6992, 1.9289)
# u3 <- c(-0.7019, -0.4083, -0.9936, -0.5439, -0.3921)
# y <- c(u1, u2, u3)
# g <- as.factor(c(rep(1, 5), rep(2, 5), rep(3, 5)))
# df0 <- data.frame(x = g, y = y)
# set.seed(2627)
# kSamples::ad.test(u1, u2, u3, method = "exact", dist = FALSE, Nsim = 1000)






