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

summarize <- function(data, formula)
{
    df0 <- stats::model.frame(formula = formula, data = data)
    colnames(df0) <- c("y", "x")
    df0 <- df0[!is.na(df0$y), , drop = FALSE]

    n_fct_lvl <- length(unique(df0$x))

    desc_mat <- with(
        data = df0,
        expr = vapply(
            X = c("sum", "length", "mean", "sd", "median", "min", "max"),
            FUN = function(fns) tapply(y, x, fns),
            FUN.VALUE = numeric(n_fct_lvl)
        )
    )

    return(desc_mat)
}

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Formula ====
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
formula_response <- function(formula)
{
    f_terms <- stats::terms(stats::as.formula(formula))
    response <- as.character(attr(f_terms, "variables"))[2]
    return(response)
}

formula_factor <- function(formula)
{
    f_terms <- stats::terms(stats::as.formula(formula))
    fcts <- as.character(attr(f_terms, "variables"))[-(1:2)]
    return(fcts)
}

formula_term <- function(formula)
{
    f_terms <- stats::terms(stats::as.formula(formula))
    ret <- as.character(attr(f_terms, "term.labels"))
    return(ret)
}

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Data types ====
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#' Check if the character could be numeric
#' Check if the character could be numeric
#'
#' @param x A character vector
#' @param include_L_sign (Default: TRUE). If TRUE, for example, "5L" will be considered as a number.
#' @param exclude_infinite (Default: TRUE). If TRUE, for example, "Inf" will NOT be considered as a number.
#' @param exclude_boolean (Default: TRUE). If TRUE, for example, "TRUE" will NOT be considered as a number.
#'
#' @return A boolean vector
#' @export
#'
#' @examples
#' x <- c("1a", ".7", "46", "2.3.3", "1.2r", "1.2", "1e4", "a1", "5L", "-.22", -Inf, NA, NaN)
#' string_is_number(x)
#' #>    1a    .7    46 2.3.3  1.2r   1.2   1e4    a1    5L  -.22  -Inf    NA   NaN
#' #> FALSE  TRUE  TRUE FALSE FALSE  TRUE  TRUE FALSE  TRUE  TRUE  FALSE FALSE FALSE
string_is_number <- function(
        x,
        include_L_sign = TRUE,
        exclude_infinite = TRUE,
        exclude_boolean = TRUE
){
    if (!is.character(x)) warning("Input should be a string vector. Coerced to character.")

    `__x__` <- x
    x <- as.character(x)

    if (include_L_sign) {
        x <- vapply(
            x,
            function(`_`) {
                if (grepl("^\\d+L$", `_`)) {
                    return(gsub("L", "", `_`))
                } else {
                    return(`_`)
                }
            },
            FUN.VALUE = character(1)
        )
        x <- unname(x)
    }

    if (exclude_infinite)
        x[is.infinite(x) | x %in% c("Inf", "-Inf")] <- NA_character_

    if (exclude_boolean) {
        x[isTRUE(x) | isFALSE(x) | x %in% c("TRUE", "FALSE")] <- NA_character_
    } else {
        x[x %in% c("TRUE", "FALSE")] <- "9"
    }

    ret <- suppressWarnings({
        !is.na(as.numeric(x))
    })
    ret <- stats::setNames(ret, `__x__`)

    return(ret)
    # ret <- vapply(
    #     X = x,
    #     FUN = function(`_x_`) {
    #         out <- try(eval(parse(text = `_x_`)), silent = TRUE)
    #         out <- !inherits(out, "try-error")
    #         return(out)
    #     },
    #     FUN.VALUE = logical(1)
    # )
    # ret[is.na(names(ret))] <- FALSE
    # ret[names(ret) %in% c("c", "C", "t", "F", "T")] <- FALSE
    # names(ret)[is.na(names(ret))] <- "NA"
    # ret[names(ret) == "NaN"] <- FALSE

    test1 <- c(
        "12", "-45.27", "0.19", "-.68", "000", "3e5", "60L",
        "NA", "NaN", "Inf", "-Inf", "TRUE", "FALSE",
        "a", "b", "c", "d", "t", "qwe", "q", "n",
        "!", "|>", "$", "#", " ", "", "&", "|", "+", "-", "*", "\\"
    )
    string_is_number(test1)

    test2 <- c(12, -45.27, 0.19, .68, 000, 3e5, 60L,
               NA, NaN, Inf, -Inf, TRUE, FALSE)
    setNames(string_is_number(test2), test2)
}


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Data structure ====
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#' Check vector
#' Check if an object is a vector.
#'
#' @param x An object
#'
#' @return Boolean value
#' @export
#'
#' @examples
#' is_vector(c(1, 2, 3))
#' #> TRUE
#'
#' is_vector(list(a = c(1, 2, 3), b = c("a", "b", "c")))
#' #> FALSE
is_vector <- function(x) { base::is.null(base::dim(x)) & base::is.atomic(x) }

is_list <- function(x) { base::is.list(x) & !base::is.atomic(x) & "list" %in% base::class(x) }

is_dataframe <- function(x) { base::length(base::dim(x)) > 1 & base::is.data.frame(x) }


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Groups homogeneity ====
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

## Check if the sample sizes of each groups are not equal
is_unbalance <- function(data, formula)
{
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


## Check if the data is tied (values are identical)
is_tied <- function(data, formula)
{
    if (is.data.frame(data)) {
        y_name <- as.character(formula)[2]
        x_name <- as.character(formula)[3]
        df0 <- data.frame(
            y = data[[y_name]],
            x = data[[x_name]]
        )
        df0 <- df0[!is.na(df0$y), ]
        res <- any(unlist(tapply(df0$y, df0$x, duplicated)))
    }

    if (is.null(dim(data))){
        res <- ( length(unique(data)) < length(data) )
    }

    return(res)
}


dataframe_to_list <- function(data, formula)
{
    if (!is.data.frame(data)) stop("Input `data` should be a dataframe")

    df0 <- stats::model.frame(formula = formula, data = data)
    colnames(df0) <- c("y", "x")

    group_names <- unique(df0[["x"]])

    lst <- list()
    for (g in group_names)
    {
        lst[[g]] <- with(
            data = df0,
            expr = subset(
                subset = (x == g),
                select = y,
                drop = TRUE
            )
        )
    }

    return(lst)
}


list_to_dataframe <- function(data, formula = NULL)
{
    lst <- data
    max_n <- base::max(base::sapply(lst, length))

    if (!is_list(lst)) stop("Input `data` should be a list")
    if (base::is.null(base::names(lst))) base::names(lst) <- base::seq_along(lst)

    if (base::is.null(formula))
    {
        x_name <- "groups"
        y_name <- "values"
    } else {
        x_name <- base::as.character(formula)[3]
        y_name <- base::as.character(formula)[2]
    }

    df0 <- base::list2DF(lst)

    df0 <- stats::reshape(
        data = df0,
        direction = "long",
        timevar = x_name,
        times = colnames(df0),
        varying = colnames(df0),
        v.names = y_name,
        idvar = "|#|@| id |&|*|"
    )

    ret <- df0[, colnames(df0) != "|#|@| id |&|*|"]

    return(ret)
}


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Equal to `outer(X = x1, Y = x1, FUN = function(x1, x2) paste(x1, x2, sep = " "))`
## applying the `paste` function to the two identical matrices
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
outer2 <- function(x, FUN = "paste", drop = TRUE)
{
    res <- base::outer(x, x, FUN)
    if (drop) res <- res[base::upper.tri(res)]
    return(res)
}


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Estimate the compact letter display (CLD) position to show on the plot
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#' Estimate the suitable compact letter y-position used in the ggplot boxplot
#' @param x The maximum values of each groups
#' @return numeric values (y-axis position)
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_point geom_text
#' @export
#' @examples
#' library(ggplot2)
#' data("iris")
#' res <- oneway_test(iris, Sepal.Length ~ Species)
#' res <- res$result
#' res$y_pos <- estimate_cld_pos(res$MAX)
#' print(res)
#' #>                 GROUP  N   AVG        SD MED MIN MAX CLD y_pos
#' #> virginica   virginica 50 6.588 0.6358796 6.5 4.9 7.9   a 9.286
#' #> versicolor versicolor 50 5.936 0.5161711 5.9 4.9 7.0   b 8.386
#' #> setosa         setosa 50 5.006 0.3524897 5.0 4.3 5.8   c 7.186
#'
#' ggplot(iris, aes(Species, Sepal.Length, color = Species)) +
#'     geom_boxplot() +
#'     geom_point() +
#'     geom_text(
#'         inherit.aes = FALSE,
#'         data = res,
#'         mapping = aes(GROUP, y_pos, label = CLD)
#'     )
estimate_cld_pos <- function(x)
{
    x <- stats::na.omit(x)
    MAX <- max(x)
    letter_pos <- x + ((ceiling(max(MAX) * 1.12) - max(x)) * 0.57)
    return(letter_pos)
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


