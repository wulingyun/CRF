\name{infer}
\alias{infer}
\alias{infer.exact}
\alias{infer.chain}
\alias{infer.tree}
\alias{infer.conditional}
\alias{infer.cutset}
\alias{infer.sample}
\alias{infer.lbp}
\alias{infer.trbp}

\title{
Inference methods
}

\description{
Computing the partition function and marginal probabilities
}

\usage{
infer.exact(crf)
infer.chain(crf)
infer.tree(crf)
infer.conditional(crf, clamped, infer.method, ...)
infer.cutset(crf, cutset, engine = "default")
infer.sample(crf, sample.method, ...)
infer.lbp(crf, max.iter = 10000, cutoff = 1e-4, verbose = 0)
infer.trbp(crf, max.iter = 10000, cutoff = 1e-4, verbose = 0)
}

\arguments{
  \item{crf}{}
  \item{clamped}{}
  \item{cutset}{}
  \item{infer.method}{}
  \item{sample.method}{}
  \item{engine}{}
  \item{max.iter}{}
  \item{cutoff}{}
  \item{verbose}{}
  \item{...}{}
}

\details{
\itemize{
  \item \code{\link{infer.exact}}           Exact inference for small graphs with brute-force counting
  \item \code{\link{infer.chain}}           Exact inference for chain-structured graphs with the forward-backward algorithm
  \item \code{\link{infer.tree}}            Exact inference for tree- and forest-structured graphs with sum-product belief propagation
  \item \code{\link{infer.conditional}}     Conditional inference (takes another inference method as input)
  \item \code{\link{infer.cutset}}          Exact inference for graphs with a small cutset using cutset conditioning
  \item \code{\link{infer.sample}}          Approximate inference using sampling (takes a sampling method as input)
  \item \code{\link{infer.lbp}}             Approximate inference using sum-product loopy belief propagation
  \item \code{\link{infer.trbp}}            Approximate inference using sum-product tree-reweighted belief propagation
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
\code{\link{make.crf}}, \code{\link{decode}}, \code{\link{sample}}
}

\examples{
library(CRF)
data(Small)
i <- infer.exact(small.crf)
i <- infer.chain(small.crf)
i <- infer.tree(small.crf)
i <- infer.conditional(small.crf, c(0,1,0,0), infer.exact)
i <- infer.cutset(small.crf, c(2))
i <- infer.sample(small.crf, sample.exact, 10000)
i <- infer.lbp(small.crf)
i <- infer.trbp(small.crf)
}

\keyword{ infer }
