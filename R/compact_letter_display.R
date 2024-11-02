#' Compact Letter Display (CLD)
#'
#' Represent significance statements resulting from all-pairwise comparisons.
#'
#' @param groups A character vector with the names of each groups (i.e. independent variables).
#' @param means A numeric vector with the mean values of each groups (i.e. response variable).
#' @param comparisons A character vector with the comparisons annotation, yielded from post-hoc test.
#' @param pval A numeric vector with the p-values corresponding to the `comparisons`.
#' @param alpha Default: 0.05; a numeric scalar.
#' @param comparison_symbol Default: " |vs| "
#' @param display_symbols Default: c("a", "b", "c", ...)
#' @param descending Default: TRUE
#'
#' @importFrom utils data
#' @importFrom stats TukeyHSD
#'
#' @return A character named-vector
#'
#' @export
#' @author Joon-Keat Lai
#'
#' @examples
#' utils::data("iris", package = "datasets")
#' stats_desc <- tapply(iris$Sepal.Length, iris$Species, "mean")
#' stats_res <- stats::TukeyHSD(aov(Sepal.Length ~ Species, iris))
#' stats_res <- as.data.frame(stats_res$Species)
#' compact_letter_display(
#'     groups = names(stats_desc),
#'     means = stats_desc,
#'     comparisons = row.names(stats_res),
#'     pval = stats_res$`p adj`,
#'     comparison_symbol = "-"
#' )
#' #> setosa versicolor  virginica
#' #>   "c"      "b"        "a"
#' @references
#' Piepho, H.P., 2004.
#' An Algorithm for a Letter-Based Representation of All-Pairwise Comparisons.
#' Journal of Computational and Graphical Statistics 13, 456–466.
#' https://doi.org/10.1198/1061860043515

compact_letter_display <- function(
        groups,
        means,
        comparisons,
        pval,
        alpha = 0.05,
        comparison_symbol = " |vs| ",
        display_symbols = base::letters,
        descending = TRUE
){
    if (!isTRUE(descending) & !isFALSE(descending)) descending <- TRUE

    if (any(grepl(comparison_symbol, groups)))
        stop("Comparison symbol should not exists in group names")

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ## The p-value of each comparisons ====
    ## the order does not matter
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    comparisons <- strsplit(comparisons, comparison_symbol, fixed = TRUE)
    comparisons <- do.call(rbind.data.frame, comparisons)
    x1 <- comparisons[, 1, drop = TRUE]
    x2 <- comparisons[, 2, drop = TRUE]


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ## Order by the mean value of each groups ====
    ## the group names will be matched to the comparisons strings to get the corresponding p-value
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    sorted_groups <- stats::setNames(means, groups)
    sorted_groups <- names(sort(sorted_groups, decreasing = descending))


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ## Generate a NULL matrix ====
    ## the row and col names are sorted by the mean-values of the groups
    ## the NULL matrix is filled with NAs
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    null_mat <- matrix(
        nrow = length(sorted_groups),
        ncol = length(sorted_groups),
        dimnames = list(sorted_groups, sorted_groups)
    )


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ## Fill in p-values ====
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    pval_mat <- null_mat
    for (i in seq_along(x1)){
        pval_mat[x1[i], x2[i]] <- pval[i]
    }


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ## Insertion step ====
    ## Non-significant comparisons were marked as TRUE,
    ## So later, non-significant comparing-pairs will be inserted with same letters.
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    bool_mat <- t(pval_mat >= alpha)
    bool_mat[is.na(bool_mat)] <- TRUE


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ## Absorption preparation ====
    ## NOT performed yet, just memorizing which columns have to be discarded
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    is_redundant <- c()
    for (n in 1:(ncol(bool_mat))){
        # m1: the nth column of the bool_mat (i.e., mat_n)
        # m2: the columns before the nth column of the bool_mat (i.e., mat_1:(n-1))
        m1 <- bool_mat[, n]
        m2 <- bool_mat[, 1:(n-1), drop = FALSE]
        is_redundant <- append(
            is_redundant,
            any(apply(m2, 2, function(x) identical(m1, x)))
        )
        # Note: When n = 1, the first column of the `m2` is always identical with the `mat-1`,
        # so, coerce the first element of the `is_redundant` vector to `FALSE` after the loop end.
    }
    is_redundant[1] <- FALSE

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ## Remove duplicate comparison-pairs
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    bool_mat[upper.tri(bool_mat)] <- FALSE


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ## Absorption step ====
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    bool_mat <- bool_mat[, !is_redundant]


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ## Substitute with letters ====
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    letter_mat <- bool_mat
    for (i in 1:ncol(bool_mat)){
        insert_letter <- bool_mat[, i, drop = TRUE]
        letter_mat[, i] <- ifelse(
            test = insert_letter,
            yes = display_symbols[i],
            no = ""
        )
    }


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ## Output
    ## The matrix will be reduced to a named-vector
    ## after row-wise collapsing the letter-columns
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    res <- apply(letter_mat, 1, function(x) paste(x, collapse = ""))
    res <- res[groups]

    return(res)
}


cld <- compact_letter_display











