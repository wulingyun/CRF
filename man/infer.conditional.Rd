\name{infer.conditional}
\alias{infer.conditional}
\title{Conditional inference method}
\usage{
infer.conditional(crf, clamped, infer.method, ...)
}
\arguments{
  \item{crf}{}

  \item{clamped}{}

  \item{infer.method}{}

  \item{...}{}
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
Conditional inference (takes another inference method as
input)
}
\examples{
library(CRF)
data(Small)
i <- infer.conditional(small.crf, c(0,1,0,0), infer.exact)
}

