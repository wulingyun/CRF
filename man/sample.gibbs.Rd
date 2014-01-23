\name{sample.gibbs}
\alias{sample.gibbs}
\title{Sampling method using single-site Gibbs sampler}
\usage{
sample.gibbs(crf, size, burn.in = 1000, start = apply(crf$node.pot, 1,
  which.max))
}
\arguments{
  \item{crf}{}

  \item{size}{}

  \item{burn.in}{}

  \item{start}{}
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
Approximate sampling using a single-site Gibbs sampler
}
\examples{
library(CRF)
data(Small)
s <- sample.gibbs(small.crf, 100)
}

