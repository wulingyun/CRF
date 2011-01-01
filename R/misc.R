make.crf <- function(adj.matrix, nstates)
{
	data <- list()
	if (!is.matrix(adj.matrix) || dim(adj.matrix)[1] != dim(adj.matrix)[2])
		stop("'adj.matrix' should be a square matrix")
	data$n.nodes <- dim(adj.matrix)[1]

	e <- which((adj.matrix + t(adj.matrix)) != 0, arr.ind = TRUE)
	data$edges <- matrix(e[e[,1] < e[,2],], ncol=2)
	data$n.edges <- nrow(data$edges)

	data$neibor.edges <- list()
	data$neibor.nodes <- list()
	for (i in 1:data$n.nodes)
	{
		data$neibor.edges[[i]] <- numeric(0)
		data$neibor.nodes[[i]] <- numeric(0)
	}
	for (i in 1:data$n.edges)
	{
		n1 <- data$edges[i, 1]
		n2 <- data$edges[i, 2]
		data$neibor.edges[[n1]] <- c(data$neibor.edges[[n1]], i)
		data$neibor.edges[[n2]] <- c(data$neibor.edges[[n2]], i)
		data$neibor.nodes[[n1]] <- c(data$neibor.nodes[[n1]], n2)
		data$neibor.nodes[[n2]] <- c(data$neibor.nodes[[n2]], n1)
	}
	for (i in 1:data$n.nodes)
	{
		data$neibor.edges[[i]] <- sort(data$neibor.edges[[i]])
		data$neibor.nodes[[i]] <- sort(data$neibor.nodes[[i]])
	}

	data$n.states <- rep(nstates, length.out = data$n.nodes)
	data$max.state <- max(nstates)

	data$node.pot <- array(1, dim=c(data$n.nodes, data$max.state))
	data$edge.pot <- array(1, dim=c(data$max.state, data$max.state, data$n.edges))

	data
}
