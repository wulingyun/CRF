\name{decode.chain}
\alias{decode.chain}
\title{Decoding method for chain-structured graphs}
\usage{
decode.chain(crf)
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
Exact decoding for chain-structured graphs with the Viterbi
algorithm.
}
\examples{
library(CRF)
data(Small)
d <- decode.chain(small.crf)
}

