\docType{package}
\name{CRF-package}
\alias{CRF}
\alias{CRF-package}
\title{CRF - Conditional Random Fields}
\description{
Library of Conditional Random Fields model
}
\details{
CRF is R package for various computational tasks of
conditional random fields as well as other probabilistic
undirected graphical models of discrete data with pairwise
and unary potentials. The decoding/inference/sampling tasks
are implemented for general discrete undirected graphical
models with pairwise potentials. The training task is less
general, focusing on conditional random fields with
log-linear potentials and a fixed structure. The code is
written entirely in R and C++. The initial version is
ported from UGM written by Mark Schmidt.

Decoding: Computing the most likely configuration \itemize{
\item \code{\link{decode.exact}} Exact decoding for small
graphs with brute-force search \item
\code{\link{decode.chain}} Exact decoding for
chain-structured graphs with the Viterbi algorithm \item
\code{\link{decode.tree}} Exact decoding for tree- and
forest-structured graphs with max-product belief
propagation \item \code{\link{decode.conditional}}
Conditional decoding (takes another decoding method as
input) \item \code{\link{decode.cutset}} Exact decoding for
graphs with a small cutset using cutset conditioning \item
\code{\link{decode.junction}} Exact decoding for
low-treewidth graphs using junction trees \item
\code{\link{decode.sample}} Approximate decoding using
sampling (takes a sampling method as input) \item
\code{\link{decode.marginal}} Approximate decoding using
inference (takes an inference method as input) \item
\code{\link{decode.lbp}} Approximate decoding using
max-product loopy belief propagation \item
\code{\link{decode.trbp}} Approximate decoding using
max-product tree-reweighted belief propagtion \item
\code{\link{decode.greedy}} Approximate decoding with
greedy algorithm \item \code{\link{decode.icm}} Approximate
decoding with the iterated conditional modes algorithm
\item \code{\link{decode.block}} Approximate decoding with
the block iterated conditional modes algorithm \item
\code{\link{decode.ilp}} Exact decoding with an integer
linear programming formulation and approximate using LP
relaxation }

Inference: Computing the partition function and marginal
probabilities \itemize{ \item \code{\link{infer.exact}}
Exact inference for small graphs with brute-force counting
\item \code{\link{infer.chain}} Exact inference for
chain-structured graphs with the forward-backward algorithm
\item \code{\link{infer.tree}} Exact inference for tree-
and forest-structured graphs with sum-product belief
propagation \item \code{\link{infer.conditional}}
Conditional inference (takes another inference method as
input) \item \code{\link{infer.cutset}} Exact inference for
graphs with a small cutset using cutset conditioning \item
\code{\link{infer.junction}} Exact decoding for
low-treewidth graphs using junction trees \item
\code{\link{infer.sample}} Approximate inference using
sampling (takes a sampling method as input) \item
\code{\link{infer.lbp}} Approximate inference using
sum-product loopy belief propagation \item
\code{\link{infer.trbp}} Approximate inference using
sum-product tree-reweighted belief propagation }

Sampling: Generating samples from the distribution
\itemize{ \item \code{\link{sample.exact}} Exact sampling
for small graphs with brute-force inverse cumulative
distribution \item \code{\link{sample.chain}} Exact
sampling for chain-structured graphs with the
forward-filter backward-sample algorithm \item
\code{\link{sample.tree}} Exact sampling for tree- and
forest-structured graphs with sum-product belief
propagation and backward-sampling \item
\code{\link{sample.conditional}} Conditional sampling
(takes another sampling method as input) \item
\code{\link{sample.cutset}} Exact sampling for graphs with
a small cutset using cutset conditioning \item
\code{\link{sample.junction}} Exact sampling for
low-treewidth graphs using junction trees \item
\code{\link{sample.gibbs}} Approximate sampling using a
single-site Gibbs sampler }

Training: Given data, computing the most likely estimates
of the parameters \itemize{ \item \code{\link{train.crf}}
Train CRF model \item \code{\link{train.mrf}} Train MRF
model }

Tools: Tools for building and manipulating CRF data
\itemize{ \item \code{\link{make.crf}} Generate CRF from
the adjacent matrix \item \code{\link{make.features}} Make
the data structure of CRF features \item
\code{\link{make.par}} Make the data structure of CRF
parameters \item \code{\link{duplicate.crf}} Duplicate an
existing CRF \item \code{\link{clamp.crf}} Generate clamped
CRF by fixing the states of some nodes \item
\code{\link{clamp.reset}} Reset clamped CRF by changing the
states of clamped nodes \item \code{\link{sub.crf}}
Generate sub CRF by selecting some nodes \item
\code{\link{mrf.update}} Update node and edge potentials of
MRF model \item \code{\link{crf.update}} Update node and
edge potentials of CRF model }
}
\examples{
library(CRF)
data(Small)
decode.exact(Small$crf)
infer.exact(Small$crf)
sample.exact(Small$crf, 100)
}
\author{
Ling-Yun Wu \email{wulingyun@gmail.com}
}
\references{
J. Lafferty, A. McCallum, and F. Pereira. Conditional
random fields: Probabilistic models for segmenting and
labeling sequence data. In \emph{the proceedings of
International Conference on Machine Learning (ICML)}, pp.
282-289, 2001.

Mark Schmidt. UGM: Matlab code for undirected graphical
models.
\url{http://www.di.ens.fr/~mschmidt/Software/UGM.html}
}
\keyword{package}

