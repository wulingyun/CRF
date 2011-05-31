\name{decode}
\alias{decode}
\alias{decode.exact}
\alias{decode.chain}
\alias{decode.tree}
\alias{decode.conditional}
\alias{decode.sample}
\alias{decode.marginal}
\alias{decode.lbp}

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
decode.sample(crf, sample.method, ...)
decode.marginal(crf, infer.method, ...)
decode.lbp(crf, max.iter=10000, cutoff=1e-4, verbose=0)
}

\arguments{
  \item{crf}{
%%     ~~Describe \code{crf} here~~
}
}

\details{
\itemize{
  \item \code{\link{decode.exact}}          Exact decoding for small graphs with brute-force search
  \item \code{\link{decode.chain}}          Exact decoding for chain-structured graphs with the Viterbi algorithm
  \item \code{\link{decode.tree}}           Exact decoding for tree- and forest-structured graphs with max-product belief propagation
  \item \code{\link{decode.conditional}}    Conditional decoding (takes another decoding method as input)
  \item \code{\link{decode.sample}}         Approximate decoding using sampling (takes a sampling method as input)
  \item \code{\link{decode.marginal}}       Approximate decoding using inference (takes an inference method as input)
  \item \code{\link{decode.lbp}}            Approximate decoding using max-product loopy belief propagation
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
decode.exact(small.crf)
decode.chain(small.crf)
decode.tree(small.crf)
decode.conditional(small.crf, c(0,1,0,0), decode.exact)
decode.sample(small.crf, sample.exact, 10000)
decode.marginal(small.crf, infer.exact)
decode.lbp(small.crf)
}

\keyword{ decode }
