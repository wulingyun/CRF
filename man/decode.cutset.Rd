\name{decode.cutset}
\alias{decode.cutset}
\title{Decoding method for graphs with a small cutset}
\usage{
decode.cutset(crf, cutset, engine = "default", start = apply(crf$node.pot,
  1, which.max))
}
\arguments{
  \item{crf}{}

  \item{cutset}{}

  \item{engine}{The underlying engine for cutset decoding,
  possible values are "default", "none", "exact", "chain",
  and "tree".}

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
Exact decoding for graphs with a small cutset using cutset
conditioning
}
\examples{
library(CRF)
data(Small)
d <- decode.cutset(small.crf, c(2))
}

