\name{crf.update}
\alias{crf.update}
\title{Update CRF potentials}
\usage{
crf.update(crf, node.fea = NaN, edge.fea = NaN, node.ext = NaN,
  edge.ext = NaN)
}
\arguments{
  \item{crf}{The CRF}

  \item{node.fea}{The node features matrix with dimension
  \code{(n.nf, n.nodes)}}

  \item{edge.fea}{The edge features matrix with dimension
  \code{(n.ef, n.edges)}}

  \item{node.ext}{The extended information of node
  features}

  \item{edge.ext}{The extended information of edge
  features}
}
\value{
This function will directly modify the CRF and return the
same CRF.
}
\description{
Update node and edge potentials of CRF model
}
\details{
This function updates \code{node.pot} and \code{edge.pot}
of CRF model by using the current values of parameters and
features.

There are two ways to model the relationship between
parameters and features. The first one exploits the special
structure of features to reduce the memory usage. However
it may not suitable for all circumstances. The other one is
more straighforward by explicitly specifying the
coefficients of each parameter to calculate the potentials,
and may use much more memory. Two approaches can be used
together.

The first way uses the objects \code{node.par} and
\code{edge.par} to define the structure of features and
provides the feature information in variables
\code{node.fea} and \code{edge.fea}. The second way
directly provides the feature information in variables
\code{node.ext} and \code{edge.ext} without any prior
assumption on feature structure. \code{node.ext} is a list
and each element has the same structure as \code{node.pot}.
\code{edge.ext} is a list and each element has the same
structure as \code{edge.pot}.

In detail, the node potential is updated as follows:

\deqn{ node.pot[n,i] = \sum_{f} par[node.par[n,i,f]] *
node.fea[f,n] + \sum_{k} par[k] * node.ext[[k]][n,i] }

and the edge potential is updated as follows:

\deqn{ edge.pot[[e]][i,j] = \sum_{f}
par[edge.par[[e]][i,j,f]] * edge.fea[f,e] + \sum_{k} par[k]
* edge.ext[[k]][[e]][i,j] }
}
\seealso{
\code{\link{crf.nll}}, \code{\link{train.crf}}
}

