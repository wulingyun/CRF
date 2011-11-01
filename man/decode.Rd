\name{decode}
\alias{decode}
\alias{decode.exact}
\alias{decode.chain}
\alias{decode.tree}
\alias{decode.conditional}
\alias{decode.cutset}
\alias{decode.sample}
\alias{decode.marginal}
\alias{decode.lbp}
\alias{decode.trbp}
\alias{decode.greedy}
\alias{decode.icm}
\alias{decode.block}

\title{
Decoding methods
}

\description{
Computing the most likely configuration
}

\usage{
decode.exact(crf)
decode.chain(crf)
decode.tree(crf)
decode.conditional(crf, clamped, decode.method, ...)
decode.cutset(crf, cutset, engine = "default",
              start = apply(crf$node.pot, 1, which.max))
decode.sample(crf, sample.method, ...)
decode.marginal(crf, infer.method, ...)
decode.lbp(crf, max.iter = 10000, cutoff = 1e-4, verbose = 0)
decode.trbp(crf, max.iter = 10000, cutoff = 1e-4, verbose = 0)
decode.greedy(crf, restart = 0, start = apply(crf$node.pot, 1, which.max))
decode.icm(crf, restart = 0, start = apply(crf$node.pot, 1, which.max))
decode.block(crf, blocks, decode.method = decode.tree, restart = 0,
             start = apply(crf$node.pot, 1, which.max), ...)
}

\arguments{
  \item{crf}{}
  \item{clamped}{}
  \item{cutset}{}
  \item{decode.method}{}
  \item{infer.method}{}
  \item{sample.method}{}
  \item{engine}{}
  \item{start}{}
  \item{max.iter}{}
  \item{cutoff}{}
  \item{verbose}{}
  \item{restart}{}
  \item{blocks}{}
  \item{...}{}
}

\details{
\itemize{
  \item \code{\link{decode.exact}}          Exact decoding for small graphs with brute-force search
  \item \code{\link{decode.chain}}          Exact decoding for chain-structured graphs with the Viterbi algorithm
  \item \code{\link{decode.tree}}           Exact decoding for tree- and forest-structured graphs with max-product belief propagation
  \item \code{\link{decode.conditional}}    Conditional decoding (takes another decoding method as input)
  \item \code{\link{decode.cutset}}         Exact decoding for graphs with a small cutset using cutset conditioning
  \item \code{\link{decode.sample}}         Approximate decoding using sampling (takes a sampling method as input)
  \item \code{\link{decode.marginal}}       Approximate decoding using inference (takes an inference method as input)
  \item \code{\link{decode.lbp}}            Approximate decoding using max-product loopy belief propagation
  \item \code{\link{decode.trbp}}           Approximate decoding using max-product tree-reweighted belief propagtion
  \item \code{\link{decode.greedy}}         Approximate decoding with greedy algorithm
  \item \code{\link{decode.icm}}            Approximate decoding with the iterated conditional modes algorithm
  \item \code{\link{decode.block}}          Approximate decoding with the block iterated conditional modes algorithm
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
\code{\link{make.crf}}, \code{\link{infer}}, \code{\link{sample}}
}

\examples{
library(CRF)
data(Small)
d <- decode.exact(small.crf)
d <- decode.chain(small.crf)
d <- decode.tree(small.crf)
d <- decode.conditional(small.crf, c(0,1,0,0), decode.exact)
d <- decode.cutset(small.crf, c(2))
d <- decode.sample(small.crf, sample.exact, 10000)
d <- decode.marginal(small.crf, infer.exact)
d <- decode.lbp(small.crf)
d <- decode.trbp(small.crf)
d <- decode.greedy(small.crf)
d <- decode.icm(small.crf)
d <- decode.block(small.crf, list(c(1,3), c(2,4)))
}

\keyword{ decode }
