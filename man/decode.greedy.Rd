\name{decode.greedy}
\alias{decode.greedy}
\title{Decoding method using greedy algorithm}
\usage{
decode.greedy(crf, restart = 0, start = apply(crf$node.pot, 1, which.max))
}
\arguments{
  \item{crf}{}

  \item{restart}{}

  \item{start}{}
}
\value{
This function will return the most likely configuration,
which is a vector of length \code{crf$n.nodes}.
}
\description{
Computing the most likely configuration for CRF
}
\details{
Approximate decoding with greedy algorithm
}
\examples{
library(CRF)
data(Small)
d <- decode.greedy(small.crf)
}

