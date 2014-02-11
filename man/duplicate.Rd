\name{duplicate}
\alias{duplicate}
\title{Duplicate CRF}
\usage{
duplicate(crf)
}
\arguments{
  \item{crf}{The existing CRF}
}
\value{
The function will return a new CRF with copied data
}
\description{
Duplicate an existing CRF
}
\details{
This function will duplicate an existing CRF. Since CRF is
implemented as an environment, normal assignment will only
copy the pointer instead of the real data. This function
will generate a new CRF and really copy all data.
}
\seealso{
\code{\link{make.crf}}
}

