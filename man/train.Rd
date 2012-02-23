\name{train}
\alias{train}
\alias{train.mrf}
\alias{train.crf}
\alias{make.features}
\alias{make.par}
\alias{mrf.update}
\alias{crf.update}
\alias{mrf.stat}
\alias{mrf.nll}
\alias{crf.nll}
\alias{gradient}

\title{
Training methods
}

\description{
Training the CRF model to estimate the parameters
}

\usage{
train.mrf(crf, instances, trace = 0)
train.crf(crf, instances, node.fea = NaN, edge.fea = NaN, node.ext = NaN, edge.ext = NaN, trace = 0)

make.features(crf, n.nf = 1, n.ef = 1)
make.par(crf, n.par = 1)

mrf.update(crf)
crf.update(crf, node.fea = NaN, edge.fea = NaN, node.ext = NaN, edge.ext = NaN)

mrf.stat(crf, instances)
mrf.nll(par, crf, instances, infer.method = infer.chain, ...)
crf.nll(par, crf, instances, node.fea = NaN, edge.fea = NaN, node.ext = NaN, edge.ext = NaN, infer.method = infer.chain, ...)
gradient(par, crf, ...)
}

\arguments{
  \item{crf}{}
  \item{instances}{}
  \item{trace}{}
  \item{node.fea}{}
  \item{edge.fea}{}
  \item{node.ext}{}
  \item{edge.ext}{}
  \item{n.nf}{}
  \item{n.ef}{}
  \item{n.par}{}
  \item{par}{}
  \item{infer.method}{}
  \item{...}{}
}

\details{
\itemize{
  \item \code{\link{train.mrf}}             Training of MRF model
  \item \code{\link{train.crf}}             Training of CRF model
  \item \code{\link{make.features}}         Make the data structure of features need for modeling and training
  \item \code{\link{make.par}}              Make the data structure of parameters need for modeling and training
  \item \code{\link{mrf.update}}            Update node.pot and edge.pot of MRF model
  \item \code{\link{crf.update}}            Update node.pot and edge.pot of CRF model
  \item \code{\link{mrf.stat}}              Calculate the sufficient statistics of MRF model
  \item \code{\link{mrf.nll}}               Calculate the negative log likelihood of MRF model
  \item \code{\link{crf.nll}}               Calculate the negative log likelihood of CRF model
  \item \code{\link{gradient}}              Calculate the gradient of negative log likelihood of CRF model
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
\code{\link{make.crf}}, \code{\link{decode}}, \code{\link{infer}}, \code{\link{sample}}
}

\examples{
library(CRF)
}

\keyword{ train }
