% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oneway_test.R
\name{oneway_test}
\alias{oneway_test}
\title{Title}
\usage{
oneway_test(
  data,
  formula,
  alpha = 0.05,
  p_adjust_method = "holm",
  use_art = TRUE
)
}
\arguments{
\item{data}{A data frame in which the variables specified in the formula will be found.}

\item{formula}{A formula specifying the model.}

\item{alpha}{Numeric value range from 0 to 1 (default: 0.05). The error tolerance.}

\item{p_adjust_method}{Character string (default: "none"). Other options: \code{stats::p.adjust.methods}.}

\item{use_art}{Whether using aligned rank transform. Default: TRUE}
}
\value{
A list, contains omnibus test and post-hoc test results.
}
\description{
Title
}
