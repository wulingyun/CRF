\name{sample}
\alias{sample}
\alias{sample.exact}
\alias{sample.chain}
\alias{sample.tree}
\alias{sample.conditional}
\alias{sample.cutset}
\alias{sample.junction}
\alias{sample.gibbs}

\title{
Sampling methods
}

\description{
Generating samples from the distribution
}

\usage{
sample.exact(crf, size)
sample.chain(crf, size)
sample.tree(crf, size)
sample.conditional(crf, size, clamped, sample.method, ...)
sample.cutset(crf, size, cutset, engine = "default")
sample.junction(crf, size)
sample.gibbs(crf, size, burn.in = 1000,
             start = apply(crf$node.pot, 1, which.max))
}

\arguments{
  \item{crf}{}
  \item{size}{}
  \item{clamped}{}
  \item{cutset}{}
  \item{sample.method}{}
  \item{engine}{}
  \item{burn.in}{}
  \item{start}{}
  \item{...}{}
}

\details{
\itemize{
  \item \code{\link{sample.exact}}          Exact sampling for small graphs with brute-force inverse cumulative distribution
  \item \code{\link{sample.chain}}          Exact sampling for chain-structured graphs with the forward-filter backward-sample algorithm
  \item \code{\link{sample.tree}}           Exact sampling for tree- and forest-structured graphs with sum-product belief propagation and backward-sampling
  \item \code{\link{sample.conditional}}    Conditional sampling (takes another sampling method as input)
  \item \code{\link{sample.cutset}}         Exact sampling for graphs with a small cutset using cutset conditioning
  \item \code{\link{sample.junction}}       Exact sampling for low-treewidth graphs using junction trees
  \item \code{\link{sample.gibbs}}          Approximate sampling using a single-site Gibbs sampler
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\references{
%% ~put references to the literature/web site here ~
}

\note{
%%  ~~further notes~~
}

\seealso{
\code{\link{make.crf}}, \code{\link{decode}}, \code{\link{infer}}
}

\examples{
library(CRF)
data(Small)
s <- sample.exact(small.crf, 100)
s <- sample.chain(small.crf, 100)
s <- sample.tree(small.crf, 100)
s <- sample.conditional(small.crf, 100, c(0,1,0,0), sample.exact)
s <- sample.cutset(small.crf, 100, c(2))
s <- sample.junction(small.crf, 100)
s <- sample.gibbs(small.crf, 100)
}

\keyword{ sample }
