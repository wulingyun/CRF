\name{CRF-package}
\alias{CRF-package}
\alias{CRF}
\docType{package}

\title{
CRF - Conditional Random Fields
}

\description{
Library to decode/infer/sample/train Conditional Random Fields
}

\details{
CRF is R package for various computational tasks of conditional random fields as well as other probabilistic undirected graphical models of discrete data with pairwise (and unary) potentials.
The decoding/inference/sampling tasks are implemented for general discrete undirected graphical models with pairwise potentials.
The training task is less general, focusing on conditional random fields with log-linear potentials and a fixed structure.
The code is written entirely in R and C++. The initial version is ported from UGM written by Mark Schmidt.

Decoding: Computing the most likely configuration
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
  \item \code{\link{decode.ilp}}            Exact decoding with an integer linear programming formulation and approximate using LP relaxation
}

Inference: Computing the partition function and marginal probabilities
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

Sampling: Generating samples from the distribution
\itemize{
  \item \code{\link{sample.exact}}          Exact sampling for small graphs with brute-force inverse cumulative distribution
  \item \code{\link{sample.chain}}          Exact sampling for chain-structured graphs with the forward-filter backward-sample algorithm
  \item \code{\link{sample.tree}}           Exact sampling for tree- and forest-structured graphs with sum-product belief propagation and backward-sampling
  \item \code{\link{sample.conditional}}    Conditional sampling (takes another sampling method as input)
  \item \code{\link{sample.cutset}}         Exact sampling for graphs with a small cutset using cutset conditioning
  \item \code{\link{sample.gibbs}}          Approximate sampling using a single-site Gibbs sampler
}

Training: Given data, computing the most likely estimates of the parameters
}

\author{
Ling-Yun Wu \email{wulingyun@gmail.com}
}

\references{
\itemize{
  \item J. Lafferty, A. McCallum, and F. Pereira. Conditional random fields: Probabilistic models for segmenting and labeling sequence data. In \emph{the proceedings of International Conference on Machine Learning (ICML)}, pp. 282-289, 2001.
  \item Mark Schmidt. UGM: A Matlab toolbox for probabilistic undirected graphical models. \url{http://www.cs.ubc.ca/~schmidtm/Software/UGM.html}
}
}

\seealso{
\code{\link{make.crf}}, \code{\link{decode}}, \code{\link{infer}}, \code{\link{sample}}, \code{\link{train}}
}

\examples{
library(CRF)
data(Small)
decode.exact(small.crf)
infer.exact(small.crf)
sample.exact(small.crf, 100)
}

\keyword{ package }
