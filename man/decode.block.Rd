\name{decode.block}
\alias{decode.block}
\title{Decoding method using block iterated conditional modes algorithm}
\usage{
decode.block(crf, blocks, decode.method = decode.tree, restart = 0,
  start = apply(crf$node.pot, 1, which.max), ...)
}
\arguments{
  \item{crf}{}

  \item{blocks}{}

  \item{decode.method}{}

  \item{restart}{}

  \item{start}{}

  \item{...}{}
}
\value{
This function will return the most likely configuration,
which is a vector of length \code{crf$n.nodes}.
}
\description{
Computing the most likely configuration for CRF
}
\details{
Approximate decoding with the block iterated conditional
modes algorithm
}
\examples{
library(CRF)
data(Small)
d <- decode.block(small.crf, list(c(1,3), c(2,4)))
}

