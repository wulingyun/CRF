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
The function will duplicate an existing CRF. Since CRF is
implemented by using environment, normal assignment will
only copy the pointer instead of the real data. The
function will really copy all data of an existing CRF to a
new CRF.
}
\seealso{
\code{\link{make.crf}}
}

