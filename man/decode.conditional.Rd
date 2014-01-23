\name{decode.conditional}
\alias{decode.conditional}
\title{Conditional decoding method}
\usage{
decode.conditional(crf, clamped, decode.method, ...)
}
\arguments{
  \item{crf}{}

  \item{clamped}{}

  \item{decode.method}{}

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
Conditional decoding (takes another decoding method as
input)
}
\examples{
library(CRF)
data(Small)
d <- decode.conditional(small.crf, c(0,1,0,0), decode.exact)
}

