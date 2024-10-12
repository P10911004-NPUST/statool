# Test if the data variances are equal (homoscedasticity)

## There are several solutions to test for the equality of variance across groups, including:
## 1. F-test
## 2. Barlett's test
## 3. Levene's test
## 4. Fligner-Killeen test


F_test <- function(data, formula){
    cat("Not yet")
}


barlett_test <- function(data, formula){
    # If the data input is a dataframe
    if (is.data.frame(data)){
        y_name <- as.character(formula)[2]
        x_name <- as.character(formula)[3]
        # res <- bartlett.test(get(y_name) ~ get(x_name), data)

        df0 <- data.frame(
            y = data[[y_name]],
            x = data[[x_name]]
        )

        df0 <- df0[!is.na(df0$x), ]

        res <- stats::bartlett.test(formula = y ~ x, data = df0)
    } else{
        cat("Data input should be a dataframe")
    }

    return(res$p.value)
}


levene_test <- function(data, formula){
    cat("Not yet")
}


fligner_test <- function(data, formula){
    # If the data input is a dataframe
    if (is.data.frame(data)){
        y_name <- as.character(formula)[2]
        x_name <- as.character(formula)[3]
        # res <- bartlett.test(get(y_name) ~ get(x_name), data)

        df0 <- data.frame(
            y = data[[y_name]],
            x = data[[x_name]]
        )

        df0 <- df0[!is.na(df0$x), ]

        res <- stats::fligner.test(formula = y ~ x, data = df0)
    } else{
        cat("Data input should be a dataframe")
    }

    return(res$p.value)
}


