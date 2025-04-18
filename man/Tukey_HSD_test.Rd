% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Tukey_HSD_test.R
\name{Tukey_HSD_test}
\alias{Tukey_HSD_test}
\title{Tukey Honestly Significant Difference (Tukey-HSD) test}
\usage{
Tukey_HSD_test(data, formula, alpha = 0.05, p_adjust_method = "none")
}
\arguments{
\item{data}{A data frame in which the variables specified in the formula will be found.}

\item{formula}{A formula specifying the model.}

\item{alpha}{Numeric value range from 0 to 1 (default: 0.05). The error tolerance.}

\item{p_adjust_method}{A character string indicating which method to use for p-value adjustment (default: "none").
For Tukey-HSD, this option should always be "none", and the calculated padj is "Tukey-adjusted p-value".
Other options: \code{stats::p.adjust.methods} (NOT recommended).}
}
\value{
A list with four vectors.
\enumerate{
\item tests: A message showing the statistical methods applied on the dataset.
\item pre_hoc: The result of pmnibus test.
\item post_hoc: includes statistics parameters for each pairwise comparisons.
\item cld: A dataframe reporting the descriptive stats and compact letter display.
}
}
\description{
Represent significance statements resulting from all-pairwise comparisons.
}
\examples{
set.seed(1)
df0 <- data.frame(
    group = as.factor(rep(c("A", "B", "C"), each = 10)),
    norm_data = c(rnorm(10, -1, 2), rnorm(10, 3, 2), rnorm(10, 0, 2)),
    skew_data = c(sqrt(rnorm(10, -7.2, 2) ^ 2), runif(10), runif(10))
)

out <- Tukey_HSD_test(df0, norm_data ~ group)
out$cld
#>   GROUP  N        AVG       SD         MED       MIN      MAX CLD
#> 1     B 10  3.4976899 2.139030  3.98374456 -1.429400 6.023562   a
#> 2     C 10 -0.2673465 1.911215  0.01843624 -3.978703 1.837955   b
#> 3     A 10 -0.7355944 1.561172 -0.48684890 -2.671257 2.190562   b
}
\references{
Howell, D.C. 2013.
Statistical methods for psychology (8th ed.). pg. 393.
Wadsworth Cengage Learning, Belmont, CA.

Tukey, John W. 1949.
Comparing Individual Means in the Analysis of Variance.
Biometrics, vol. 5, no. 2, pp. 99-114. JSTOR.
https://doi.org/10.2307/3001913.
}
\author{
Joon-Keat Lai
}
