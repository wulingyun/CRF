\name{decode.lbp}
\alias{decode.lbp}
\title{Decoding method using loopy belief propagation}
\usage{
decode.lbp(crf, max.iter = 10000, cutoff = 1e-04, verbose = 0)
}
\arguments{
  \item{crf}{}

  \item{max.iter}{}

  \item{cutoff}{}

  \item{verbose}{}
}
\value{
This function will return the most likely configuration,
which is a vector of length \code{crf$n.nodes}.
}
\description{
Computing the most likely configuration for CRF
}
\details{
Approximate decoding using max-product loopy belief
propagation
}
\examples{
library(CRF)
data(Small)
d <- decode.lbp(small.crf)
}

