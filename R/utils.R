

#' Standard deviation of the population
#'
#' @param x A numeric vector
#'
#' @return
#' standard deviation value for the population
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' sd_population(rnorm(100))
#' #> 0.8936971
sd_population <- function(x){
    # Standard deviation of population
    x <- x[stats::complete.cases(x)]
    N <- length(x)
    ret <- stats::sd(x) * sqrt( (N - 1) / N )
    return(ret)
}


#' Check vector
#'
#' Check if an object is a vector.
#'
#' @param x An object
#'
#' @return
#' Boolean value
#'
#' @export
#'
#' @examples
#' is_vector(c(1, 2, 3))
#' #> TRUE
#'
#' is_vector(list(a = c(1, 2, 3), b = c("a", "b", "c")))
#' #> FALSE
is_vector <- function(x){
    is.null(dim(x)) & is.atomic(x)
}

is_not_vector <- function(x){
    !is.null(dim(x)) & !is.atomic(x)
}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Check if sample sizes are not equal
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
is_unbalance <- function(data, formula){
    x_name <- as.character(formula)[3]
    y_name <- as.character(formula)[2]
    df0 <- data.frame(
        x = data[[x_name]],
        y = data[[y_name]]
    )
    n <- unname(with(df0, tapply(y, x, "length")))
    res <- length(unique(n)) > 1
    return(res)
}


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Check if the data is tied (all values are identical)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
is_tied <- function(data, formula){
    if (is.data.frame(data)) {
        y_name <- as.character(formula)[2]
        x_name <- as.character(formula)[3]

        df0 <- data.frame(
            y = data[[y_name]],
            x = data[[x_name]]
        )

        df0 <- df0[!is.na(df0$y), ]

        res <- any(with(df0, tapply(y, x, "stats::sd")) == 0)
    }

    if (is.null(dim(data))){
        res <- ( length(unique(data)) < length(data) )
    }

    return(res)
}



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Check data normality
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
is_normal <- function(data, formula, alpha = 0.05) {
    #### Dataframe input ====
    if (is.data.frame(data)) {
        y_name <- as.character(formula)[2]
        x_name <- as.character(formula)[3]

        df0 <- data.frame(
            y = data[[y_name]],
            x = data[[x_name]]
        )

        df0 <- df0[!is.na(df0$y), ]

        if (is_tied(df0, y ~ x)) {
            warning("This is tied-data.")
            res <- FALSE
        } else {
            aov_mod <- stats::aov(formula = y ~ x, data = df0)
            res <- stats::shapiro.test(aov_mod$residuals)$p.value > alpha
        }
    }

    #### Vector input ====
    if (is.null(dim(data))) {
        data <- stats::na.omit(data)

        if (is_tied(data)) {
            warning("This is tied-data.")
            res <- FALSE
        } else {
            res <- stats::shapiro.test(data)$p.value > alpha
        }
    }

    return(res)
}


dataframe_to_list <- function(data, formula){
    if (!is.data.frame(data)) stop("Input `data` should be a dataframe")

    df0 <- stats::model.frame(formula = formula, data = data)
    colnames(df0) <- c("y", "x")

    group_names <- unique(df0$x)

    lst <- list()
    for (g in group_names){
        lst[[g]] <- subset(df0, x == g, y, drop = TRUE)
    }
    return(lst)
}


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Equal to `outer(X = x1, Y = x1, FUN = function(x1, x2) paste(x1, x2, sep = " "))`
## applying the `paste` function to the two identical matrices
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
outer2 <- function(x, FUN = "paste", drop = TRUE){
    res <- outer(x, x, FUN)
    if (drop) res <- res[upper.tri(res)]
    return(res)
}


# qDunnett <- function (p, df, k, rho,
#                       type = c("two-sided", "one-sided"))
# {
#     type <- match.arg(type)
#     alpha <- 1 - p
#     if (type == "two-sided") {
#         alpha <- alpha/2
#     }
#     S <- matrix(rho, nrow=k, ncol=k) + (1-rho)*diag(k)
#     if (type == "two-sided") {
#         f <- function(d, df, k, S, p) {
#             mnormt::sadmvt(df=df, lower=rep(-d,k), upper=rep(d,k),
#                            mean=rep(0,k), S=S, maxpts=2000*k) - p
#         }
#     }
#     else {
#         f <- function(d, df, k, S, p) {
#             mnormt::pmt(d, S=S, df=df) - p
#         }
#     }
#     d <- uniroot(f,
#                  df = df, k = k, S = S, p=p,
#                  lower=qt(1 - alpha, df),
#                  upper=qt(1 - alpha/k, df),
#                  tol=.Machine$double.eps, maxiter=5000)$root
#     return(d)
# }


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Estimate the compact letter display (CLD) position to show on the plot
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
estimate_cld_pos <- estimate_letter_pos <- function(x){
    MAX <- max(x)
    letter_pos <- MAX + ((ceiling(max(MAX) * 1.15) - max(MAX)) * 0.43)
    return(letter_pos)
}



search_sorted <- function(x, insertion, side = "left", descending = FALSE){
    if (is_not_vector(x)) stop("Input `x` should be a vector")
    if (is_not_vector(insertion)) stop("`insertion` should be a vector")

    side <- match.arg(side, c("left", "right"))
    x <- sort(x, decreasing = descending)

    .get_ind <- function(x, insertion, side, descending){
        ind <- which(x == insertion)
        if (insertion < min(x)) return(1)
        if (insertion > max(x)) return(length(x) + 1)
        if (side == "left") return (min(ind))
        if (side == "right") return (max(ind))
        return(ind)
    }

    ret <- sapply(
        X = insertion,
        FUN = function(i) .get_ind(x = x, insertion = i, side = side, descending = descending)
    )

    return(ret)
}

























# #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# #> Plot a boxplot with compact letter display (cld)
# #> The cld should be a named vector, for example, c(level_01 = "a", level_02 = "b")
# show_boxplot <- function(
        #         data,
#         formula,
#         cld = NULL,
#         color = NULL,
#         point_size = 3,
#         triangle_size = 3,
#         letter_size = 8,
#         letter_pos_adjust = 0,
#         x_label_order = NULL
# ){
#     x_name <- as.character(formula)[3]
#     y_name <- as.character(formula)[2]
#
#     df0 <- data.frame(
#         x = data[[x_name]],
#         y = data[[y_name]]
#     )
#
#     if (is.null(x_label_order)) x_label_order <- sort(unique(as.character(df0$x)))
#     factor_missing <- any(!(x_label_order %in% unique(df0$x)))
#     if (factor_missing) warning("Factor not match")
#     df0$x <- factor(df0$x, levels = x_label_order)
#
#     if (!is.null(color)) df0$color <- data[[color]]
#
#     p1 <- ggplot(df0, aes(x, y, color = color)) +
#         theme_bw() +
#         labs(
#             x = sym(x_name),
#             y = sym(y_name)
#         ) +
#         geom_point(
#             size = point_size,
#             position = position_jitter(width = 0.1),
#             alpha = 0.5,
#             show.legend = FALSE
#         ) +
#         geom_boxplot(
#             outliers = FALSE,
#             outlier.colour = "transparent",
#             outlier.size = 0,
#             outlier.shape = NA,
#             fill = NA,
#             size = 1,
#             linewidth = 0.5
#         ) +
#         stat_summary(
#             geom = "point",
#             fun = "mean",
#             size = triangle_size,
#             color = "black",
#             shape = 17,
#             alpha = 0.7
#         ) +
#         theme(
#             text = element_text(family = "sans", face = "bold", size = 21),
#             axis.title.x.bottom = element_text(margin = ggplot2::margin(t = 9)),
#             axis.title.y.left = element_text(margin = ggplot2::margin(r = 9)),
#             legend.position = "none"
#         )
#
#     if (!is.null(cld)){
#         cld <- data.frame(
#             x = names(cld),
#             letter = unname(cld)
#         )
#
#         letter_pos <- df0 %>%
#             summarise(
#                 letter_pos = estimate_letter_pos(y),
#                 .by = x
#             ) %>%
#             left_join(cld, by = "x")
#
#         # letter_pos$letter <- unname(cld[match(names(cld), letter_pos[["x"]])])
#
#         p1 <- p1  +
#             geom_text(
#                 data = letter_pos,
#                 mapping = aes(x, letter_pos + letter_pos_adjust, label = letter),
#                 inherit.aes = FALSE,
#                 size = letter_size
#             )
#     }
#
#     return(p1)
# }


