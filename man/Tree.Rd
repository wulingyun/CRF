\docType{data}
\name{Tree}
\alias{Tree}
\title{Tree CRF example}
\format{A list containing two elements:
\itemize{
  \item \code{crf} The CRF
  \item \code{answer} A list of 4 elements:
  \itemize{
    \item \code{decode} The most likely configuration
    \item \code{node.bel} The node belief
    \item \code{edge.bel} The edge belief
    \item \code{logZ} The logarithmic value of CRF normalization factor Z
  }
}}
\usage{
data(Tree)
}
\description{
This data set gives a tree CRF example
}
\keyword{datasets}

