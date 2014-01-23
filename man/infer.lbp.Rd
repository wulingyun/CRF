\name{infer.lbp}
\alias{infer.lbp}
\title{Inference method using loopy belief propagation}
\usage{
infer.lbp(crf, max.iter = 10000, cutoff = 1e-04, verbose = 0)
}
\arguments{
  \item{crf}{}

  \item{max.iter}{}

  \item{cutoff}{}

  \item{verbose}{}
}
\value{
This function will return a list with components:
\item{node.bel}{Node belief. It is a matrix with
\code{crf$n.nodes} rows and \code{crf$max.state} columns.}
\item{edge.bel}{Edge belief. It is a list of matrices. The
size of list is \code{crf$n.edges} and the matrix \code{i}
has \code{crf$n.states[crf$edges[i,1]]} rows and
\code{crf$n.states[crf$edges[i,2]]} columns.}
\item{logZ}{The logarithmic value of CRF normalization
factor Z.}
}
\description{
Computing the partition function and marginal probabilities
}
\details{
Approximate inference using sum-product loopy belief
propagation
}
\examples{
library(CRF)
data(Small)
i <- infer.lbp(small.crf)
}

