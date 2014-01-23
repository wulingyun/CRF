\name{decode.icm}
\alias{decode.icm}
\title{Decoding method using iterated conditional modes algorithm}
\usage{
decode.icm(crf, restart = 0, start = apply(crf$node.pot, 1, which.max))
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
Approximate decoding with the iterated conditional modes
algorithm
}
\examples{
library(CRF)
data(Small)
d <- decode.icm(small.crf)
}

