\name{decode.sample}
\alias{decode.sample}
\title{Decoding method using sampling}
\usage{
decode.sample(crf, sample.method, ...)
}
\arguments{
  \item{crf}{}

  \item{sample.method}{}

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
Approximate decoding using sampling (takes a sampling
method as input)
}
\examples{
library(CRF)
data(Small)
d <- decode.sample(small.crf, sample.exact, 10000)
}

