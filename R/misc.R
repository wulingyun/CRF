make.crf <- function(adj.matrix, nstates)
{
	data <- list()
	if (!is.matrix(adj.matrix) || dim(adj.matrix)[1] != dim(adj.matrix)[2])
		stop("'adj.matrix' should be a square matrix")
	data$n.nodes <- dim(adj.matrix)[1]

	e <- which(adj.matrix != 0, arr.ind = TRUE)
	e <- matrix(c(e, e[,2], e[,1]), ncol=2)
	e <- unique(matrix(e[e[,1] < e[,2],], ncol=2))
	data$edges <- matrix(e[order(e[,1], e[,2]),], ncol=2)
	data$n.edges <- nrow(data$edges)

	data <- make.adj.info(data)

	data$n.states <- rep(nstates, length.out=data$n.nodes)
	data$max.state <- max(nstates)

	data$node.pot <- array(1, dim=c(data$n.nodes, data$max.state))
	data$edge.pot <- lapply(1:data$n.edges, function(i) array(1, dim=c(data$n.states[data$edges[i,1]], data$n.states[data$edges[i,2]])))

	class(data) <- "CRF"
	data
}

make.adj.info <- function(data)
{
	data$n.adj <- rep(0, length.out=data$n.nodes)
	data$adj.edges <- list()
	data$adj.nodes <- list()
	for (i in 1:data$n.nodes)
	{
		data$adj.edges[[i]] <- numeric(0)
		data$adj.nodes[[i]] <- numeric(0)
	}
	if (data$n.edges > 0)
	{
		for (i in 1:data$n.edges)
		{
			n1 <- data$edges[i, 1]
			n2 <- data$edges[i, 2]
			data$adj.edges[[n1]] <- c(data$adj.edges[[n1]], i)
			data$adj.edges[[n2]] <- c(data$adj.edges[[n2]], i)
			data$adj.nodes[[n1]] <- c(data$adj.nodes[[n1]], n2)
			data$adj.nodes[[n2]] <- c(data$adj.nodes[[n2]], n1)
		}
		for (i in 1:data$n.nodes)
		{
			data$n.adj[i] <- length(data$adj.edges[[i]])
			data$adj.edges[[i]] <- sort(data$adj.edges[[i]])
			data$adj.nodes[[i]] <- sort(data$adj.nodes[[i]])
		}
	}
	data
}

get.potential <- function(crf, configuration)
	.Call("Get_Potential", crf, configuration)

get.logPotential <- function(crf, configuration)
	.Call("Get_LogPotential", crf, configuration)
