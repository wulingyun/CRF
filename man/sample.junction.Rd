\name{sample.junction}
\alias{sample.junction}
\title{Sampling method for low-treewidth graphs}
\usage{
sample.junction(crf, size)
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
Exact sampling for low-treewidth graphs using junction
trees
}
\examples{
library(CRF)
data(Small)
s <- sample.junction(small.crf, 100)
}

