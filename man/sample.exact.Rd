\name{sample.exact}
\alias{sample.exact}
\title{Sampling method for small graphs}
\usage{
sample.exact(crf, size)
}
\arguments{
  \item{crf}{}

  \item{size}{}
}
\value{
This function will return a matrix with \code{size} rows
and \code{crf$n.nodes} columns, in which each row is a
sampled configuration.
}
\description{
Generating samples from the distribution
}
\details{
Exact sampling for small graphs with brute-force inverse
cumulative distribution
}
\examples{
library(CRF)
data(Small)
s <- sample.exact(small.crf, 100)
}

