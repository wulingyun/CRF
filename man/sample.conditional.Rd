\name{sample.conditional}
\alias{sample.conditional}
\title{Conditional sampling method}
\usage{
sample.conditional(crf, size, clamped, sample.method, ...)
}
\arguments{
  \item{crf}{}

  \item{size}{}

  \item{clamped}{}

  \item{sample.method}{}

  \item{...}{}
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
Conditional sampling (takes another sampling method as
input)
}
\examples{
library(CRF)
data(Small)
s <- sample.conditional(small.crf, 100, c(0,1,0,0), sample.exact)
}

