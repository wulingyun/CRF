make.crf <- function(adj.matrix, nstates)
{
	data <- new.env()
	if (!is.matrix(adj.matrix) || dim(adj.matrix)[1] != dim(adj.matrix)[2])
		stop("'adj.matrix' should be a square matrix")
	data$n.nodes <- dim(adj.matrix)[1]

	e <- which(adj.matrix != 0, arr.ind = TRUE)
	e <- matrix(c(e, e[,2], e[,1]), ncol=2)
	e <- unique(matrix(e[e[,1] < e[,2],], ncol=2))
	data$edges <- matrix(e[order(e[,1], e[,2]),], ncol=2)
	data$n.edges <- nrow(data$edges)

	.Call("Make_AdjInfo", data)

	data$n.states <- rep(nstates, length.out=data$n.nodes)
	data$max.state <- max(nstates)

	data$node.pot <- array(1, dim=c(data$n.nodes, data$max.state))
	data$edge.pot <- lapply(1:data$n.edges, function(i) array(1, dim=c(data$n.states[data$edges[i,1]], data$n.states[data$edges[i,2]])))

	class(data) <- "CRF"
	data
}

duplicate <- function(crf)
{
	data <- new.env()
	for (i in ls(envir=crf)) assign(i, get(i, envir=crf), envir=data)
	data
}

get.potential <- function(crf, configuration)
	.Call("Get_Potential", crf, configuration)

get.logPotential <- function(crf, configuration)
	.Call("Get_LogPotential", crf, configuration)
