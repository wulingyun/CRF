\name{mrf.update}
\alias{mrf.update}
\title{Update MRF potentials}
\usage{
mrf.update(crf)
}
\arguments{
  \item{crf}{The CRF}
}
\value{
This function will directly modify the CRF and return the
same CRF.
}
\description{
Update node.pot and edge.pot of MRF model
}
\details{
The function updates \code{node.pot} and \code{edge.pot} of
MRF model.
}
\seealso{
\code{\link{mrf.nll}}, \code{\link{train.mrf}}
}

