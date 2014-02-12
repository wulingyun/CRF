\name{sample.conditional}
\alias{sample.conditional}
\title{Conditional sampling method}
\usage{
sample.conditional(crf, size, clamped, sample.method, ...)
}
\arguments{
  \item{crf}{The CRF}

  \item{size}{The sample size}

  \item{clamped}{The vector of fixed values for clamped
  nodes, 0 for unfixed nodes}

  \item{sample.method}{The sampling method to solve the
  clamped CRF}

  \item{...}{The parameters for \code{sample.method}}
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
s <- sample.conditional(Small$crf, 100, c(0,1,0,0), sample.exact)
}

