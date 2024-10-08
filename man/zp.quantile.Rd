% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CalADM.R
\name{zp.quantile}
\alias{zp.quantile}
\title{Quantile Normalization Function}
\usage{
zp.quantile(x, y)
}
\arguments{
\item{x}{A numeric vector. The reference vector whose quantiles will be used for normalization.}

\item{y}{A numeric vector. The vector to be normalized based on the quantiles of `x`.}
}
\value{
A numeric vector of the same length as `y`, containing the quantile-normalized values.
}
\description{
Quantile Normalization Function
}
