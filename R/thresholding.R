


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Otsu threshold
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
otsu_threshold <- function(x) {
    raw_x <- stats::na.omit(x)
    x_sd <- 3 * stats::sd(raw_x)
    raw_x <- raw_x[ (raw_x < (raw_x + x_sd)) & (raw_x > (raw_x - x_sd)) ]
    raw_x <- sort(raw_x, decreasing = FALSE)

    x_length <- length(raw_x)

    MSE_v <- c()
    for (i in seq_along(raw_x)) {
        x1 <- raw_x[ 1 : i ]
        x2 <- raw_x[ i : x_length ]

        MSE_x1 <- sum( (x1 - mean(x1)) ^ 2 ) / length(x1)
        MSE_x2 <- sum( (x2 - mean(x2)) ^ 2 ) / length(x2)

        MSE_v[i] <- MSE_x1 + MSE_x2
    }

    MSE_v <- stats::setNames(MSE_v, as.character(raw_x))

    threshold <- MSE_v[which.min(MSE_v)]
    attr(threshold, "MSE") <- MSE_v

    return(threshold)
}
