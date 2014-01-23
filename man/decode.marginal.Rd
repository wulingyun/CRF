\name{decode.marginal}
\alias{decode.marginal}
\title{Decoding method using inference}
\usage{
decode.marginal(crf, infer.method, ...)
}
\arguments{
  \item{crf}{}

  \item{infer.method}{}

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
Approximate decoding using inference (takes an inference
method as input)
}
\examples{
library(CRF)
data(Small)
d <- decode.marginal(small.crf, infer.exact)
}

