
search_sorted <- function(x, ref, side = "left", ref_descending = FALSE, zero_index = TRUE)
{
    if (!is.null(ref_descending))
        ref <- sort(ref, decreasing = ref_descending)

    .insert_left <- function(x, ref)
    {
        seq_ref <- seq_along(ref)
        for (i in seq_ref) {
            if (ref_descending) {  # Descending order
                if (x >= max(ref)) return(1)
                if (x == min(ref)) return(length(ref))
                if (x < min(ref)) return(length(ref) + 1)
                if (x >= ref[i]) return(i)
            } else {  # Ascending order
                if (x == max(ref)) return(length(ref))
                if (x > max(ref)) return(length(ref) + 1)
                if (x < min(ref)) return(1)
                if (x <= ref[i]) return(i)
            }
        }
    }

    .insert_right <- function(x, ref)
    {
        seq_ref <- rev(seq_along(ref))
        for (i in seq_ref) {
            if (ref_descending) {  # Descending order
                if (x > max(ref)) return(1)
                if (x <= min(ref)) return(length(ref) + 1)
                if (x <= ref[i]) return(i + 1)
            } else {  # Ascending order
                if (x >= max(ref)) return(length(ref) + 1)
                if (x < min(ref)) return(1)
                if (x >= ref[i]) return(i + 1)
            }

        }
    }

    if (side == "left") ret <- vapply(x, function(v) .insert_left(v, ref), FUN.VALUE = numeric(1))
    if (side == "right") ret <- vapply(x, function(v) .insert_right(v, ref), FUN.VALUE = numeric(1))

    if (zero_index) ret <- ret - 1

    return(ret)
}





