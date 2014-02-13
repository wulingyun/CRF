\name{make.par}
\alias{make.par}
\title{Make CRF parameters}
\usage{
make.par(crf, n.par = 1)
}
\arguments{
  \item{crf}{The CRF}

  \item{n.par}{The number of parameters}
}
\value{
This function will directly modify the CRF and return the
same CRF.
}
\description{
Make the data structure of CRF parameters
}
\details{
This function makes the data structure of parameters need
for modeling and training CRF. The parameters are stored in
\code{par}, which is a numeric vector of length
\code{n.par}.
}
\seealso{
\code{\link{crf.update}}, \code{\link{make.features}},
\code{\link{make.crf}}
}

