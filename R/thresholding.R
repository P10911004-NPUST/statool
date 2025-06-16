#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Otsu threshold
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
otsu_threshold <- function(x)
{
    raw_x <- stats::na.omit(x)
    x_sd <- 3 * stats::sd(raw_x)
    raw_x <- raw_x[ (raw_x < (raw_x + x_sd)) & (raw_x > (raw_x - x_sd)) ]
    raw_x <- sort(raw_x, decreasing = FALSE)

    x_length <- length(raw_x)

    SS_v <- c()  # Sum of Squares
    for (i in seq_along(raw_x)) {
        x1 <- raw_x[ 1 : i ]
        x2 <- raw_x[ i : x_length ]

        SS_x1 <- sum( (x1 - mean(x1)) ^ 2 ) / length(x1)
        SS_x2 <- sum( (x2 - mean(x2)) ^ 2 ) / length(x2)

        SS_v[i] <- SS_x1 + SS_x2
    }

    SS_v <- stats::setNames(SS_v, as.character(raw_x))

    threshold <- SS_v[which.min(SS_v)]
    attr(threshold, "SS") <- SS_v

    return(threshold)
}
