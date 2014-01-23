\name{sample.cutset}
\alias{sample.cutset}
\title{Sampling method for graphs with a small cutset}
\usage{
sample.cutset(crf, size, cutset, engine = "default")
}
\arguments{
  \item{crf}{}

  \item{size}{}

  \item{cutset}{}

  \item{engine}{The underlying engine for cutset sampling,
  possible values are "default", "none", "exact", "chain",
  and "tree".}
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
Exact sampling for graphs with a small cutset using cutset
conditioning
}
\examples{
library(CRF)
data(Small)
s <- sample.cutset(small.crf, 100, c(2))
}

