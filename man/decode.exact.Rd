\name{decode.exact}
\alias{decode.exact}
\title{Decoding method for small graphs}
\usage{
decode.exact(crf)
}
\arguments{
  \item{crf}{}
}
\value{
This function will return the most likely configuration,
which is a vector of length \code{crf$n.nodes}.
}
\description{
Computing the most likely configuration for CRF
}
\details{
Exact decoding for small graphs with brute-force search
}
\examples{
library(CRF)
data(Small)
d <- decode.exact(small.crf)
}

