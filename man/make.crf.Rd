\name{make.crf}
\alias{make.crf}

\title{
Make CRF data structure
}

\description{
Generate CRF data structure from the adjacent matrix
}

\usage{
make.crf(adj.matrix, nstates)
}

\arguments{
  \item{adj.matrix}{The adjacent matrix of CRF network}
  \item{nstates}{The state numbers of nodes}
}

\details{
The function will generate a empty CRF data structure from a given adjacent matrix.
If the length of \code{nstates} is less than \code{n.nodes}, it will be used repeatly.
All node and edge potentials are initilized as 1.
}

\value{
The function will return a CRF data structure, which is a list with components:
\item{n.nodes}{The number of nodes.}
\item{n.edges}{The number of edges.}
\item{n.states}{The number of states for each node. It is a vector of length \code{n.nodes}.}
\item{max.state}{The maximum number of states. It is equal to \code{max(n.states)}.}
\item{edges}{The node pair of each edge. It is a matrix with 2 columns and \code{n.edges} rows. Each row denotes one edge. The node with smaller id is put in the first column.}
\item{n.adj}{The number of adjacent nodes for each node. It is a vector of length \code{n.nodes}.}
\item{adj.nodes}{The list of adjacent nodes for each node. It is a list of length \code{n.nodes} and the i-th element is a vector of length \code{n.adj[i]}.}
\item{adj.edges}{The list of adjacent edges for each node. It is similiar to \code{adj.nodes} while contains the edge ids instead of node ids.}
\item{node.pot}{The node potentials. It is a matrix with dimmension \code{(n.nodes, max.state)}. Each row \code{node.pot[i,]} denotes the node potentials of the i-th node.}
\item{edge.pot}{The edge potentials. It is a list of \code{n.edges} matrixes. Each matrix \code{edge.pot[[i]]}, with dimension \code{(n.states[edges[i,1]], n.states[edges[i,2]])}, denotes the edge potentials of the i-th edge.}
}

\seealso{
\code{\link{clamp.crf}}, \code{\link{sub.crf}}
}

\examples{
library(CRF)

nNodes <- 4
nStates <- 2

adj <- matrix(0, nrow=nNodes, ncol=nNodes)
for (i in 1:(nNodes-1))
{
	adj[i,i+1] <- 1
	adj[i+1,i] <- 1
}

crf <- make.crf(adj, nStates)

crf$node.pot[1,] <- c(1, 3)
crf$node.pot[2,] <- c(9, 1)
crf$node.pot[3,] <- c(1, 3)
crf$node.pot[4,] <- c(9, 1)

for (i in 1:crf$n.edges)
{
   crf$edge.pot[[i]][1,] <- c(2, 1)
   crf$edge.pot[[i]][2,] <- c(1, 2)
}

}

\keyword{ make.crf }
