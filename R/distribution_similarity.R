# Test whether the samples come from the same distribution

ksample_anderson_darling <- function(data = NULL, formula = NULL, ...){
    if (is.null(data) & is.null(formula)){
        data <- list(...)
    }

    df0 <- stats::model.frame(formula = formula, data = data)
    colnames(df0) <- c("y", "x")

    df0 <- df0[stats::complete.cases(df0), ]
    df0 <- df0[order(df0$y, decreasing = FALSE), ]

    desc_mat <- with(
        data = df0,
        expr = sapply(
            X = c("length", "mean"),
            FUN = function(fns) tapply(y, x, fns)
        )
    )

    group_names <- rownames(desc_mat)
    group_sizes <- desc_mat[, "length"]
    N <- nrow(df0)
    n <- length(group_names)
    if (n < 2) stop("Required at least 2 groups")

    Z <- df0$y
    Zstar <- sort(unique(Z), decreasing = FALSE)
    if (length(Zstar) < 2) stop("Needs more than one distinct observation")

    Zleft <- search_sorted(Z, Zstar, side = "left")
    Zright <- search_sorted(Z, Zstar, side = "right")

    if (N == length(Zstar)) Lj <- 1
    if (N != length(Zstar)) Lj <- Zright - Zleft

    Baj <- Zleft + (Lj / 2)

    inner <- c()
    for (g in group_names){
        ni <- group_sizes[[g]]
        xi <- subset(df0, x == g, y)[, , drop = TRUE]
        xright <- search_sorted(xi, Zstar, side = "right")
        xleft <- search_sorted(xi, Zstar, side = "left")
        fij <- xright - xleft
        Maij <- xright - (fij / 2)

        block_A <- 1 / ni
        block_B <- Lj / N
        block_C <- ( (N * Maij - ni * Baj) ^ 2 ) / ( (Baj * N) - (Baj ^ 2) - (N * Lj / 4) )
        inner <- append( inner, block_A * sum(block_B * block_C) )
    }

    AakN2 <- ((N - 1) / N) * sum(inner)

    # .AakN2 <- function(group_name){
    #     ni <- group_sizes[[group_name]]
    #     xi <- subset(df0, x == group_name, y)[, , drop = TRUE]
    #     xright <- search_sorted(xi, Zstar, side = "right")
    #     xleft <- search_sorted(xi, Zstar, side = "left")
    #     fij <- xright - xleft
    #     Maij <- xright - (fij / 2)
    #
    #     block_A <- 1 / ni
    #     block_B <- Lj / N
    #     block_C <- ( (N * Maij - ni * Baj) ^ 2 ) / ( (Baj * N) - (Baj ^ 2) - (N * Lj / 4) )
    #     inner <- block_A * sum(block_B * block_C)
    #     return(inner)
    # }
    #
    # ret <- sapply(group_names, .AakN2)


    return(AakN2)
}




# u1 <- c(1.0066, -0.9587,  0.3462, -0.2653, -1.3872)
# u2 <- c(0.1005, 0.2252, 0.4810, 0.6992, 1.9289)
# u3 <- c(-0.7019, -0.4083, -0.9936, -0.5439, -0.3921)
# y <- c(u1, u2, u3)
# g <- as.factor(c(rep(1, 5), rep(2, 5), rep(3, 5)))
# df0 <- data.frame(x = g, y = y)
# set.seed(2627)
# kSamples::ad.test(u1, u2, u3, method = "exact", dist = FALSE, Nsim = 1000)






