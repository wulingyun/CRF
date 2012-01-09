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

	adj.info <- .Call("Make_AdjInfo", data)
	data$n.adj <- adj.info$n.adj
	data$adj.nodes <- adj.info$adj.nodes
	data$adj.edges <- adj.info$adj.edges

	data$n.states <- rep(nstates, length.out=data$n.nodes)
	data$max.state <- max(nstates)

	data$node.pot <- array(1, dim=c(data$n.nodes, data$max.state))
	data$edge.pot <- lapply(1:data$n.edges, function(i) array(1, dim=c(data$n.states[data$edges[i,1]], data$n.states[data$edges[i,2]])))

	class(data) <- "CRF"
	data
}

make.features <- function(crf, n.nf = 1, n.ef = n.nf)
{
	data <- list()
	data$n.nf <- n.nf
	data$n.ef <- n.ef
	data$node.pot <- array(0, dim=c(crf$n.nodes, crf$max.state, n.nf))
	data$edge.pot <- lapply(1:crf$n.edges, function(i) array(0, dim=c(crf$n.states[crf$edges[i,1]], crf$n.states[crf$edges[i,2]], n.ef)))
	data$n.params <- 0
	data$params <- numeric(0)
	data
}

make.node.pot <- function(crf, features, instance = matrix(1, crf$n.nodes, features$n.nf))
{
	map <- features$node.pot > 0
	pot <- features$node.pot
	pot[map] <- features$params[features$node.pot[map]]
	node.pot <- exp(rowSums(pot * as.vector(apply(instance, 2, function(i) rep(i, crf$max.state))), dims=2))
}

make.edge.pot <- function(crf, features, instance = matrix(1, crf$n.edges, features$n.ef))
{
	make.edge.pot.i <- function(i)
	{
		map <- features$edge.pot[[i]] > 0
		pot <- features$edge.pot[[i]]
		pot[map] <- features$params[features$edge.pot[[i]][map]]
		pot <- exp(rowSums(pot * rep(instance[i,], each=length(crf$edge.pot[[i]])), dims=2))
	}
	edge.pot <- lapply(1:crf$n.edges, make.edge.pot.i)
}

calc.frequency <- function(v, n)
	.Call("Calc_Frequency", v, n)

mrf.suff.stat <- function(crf, features, instances)
{
	suff.stat <- numeric(features$n.params)
	n.instances <- dim(instances)[1]
	index <- cbind(rep(1:crf$n.nodes, each=n.instances), as.vector(instances), rep(1, length(instances)))
	index <- features$node.pot[index]
	suff.stat <- suff.stat + calc.frequency(index, features$n.params)
	for (i in 1:crf$n.edges)
	{
		index <- cbind(instances[,crf$edges[i,1]], instances[,crf$edges[i,2]], rep(1, n.instances))
		index <- features$edge.pot[[i]][index]
		suff.stat <- suff.stat + calc.frequency(index, features$n.params)
	}
	suff.stat
}

mrf.nll <- function(crf, features, instances, infer.method=infer.chain, ...)
{
	nll <- list()
	n.instances <- dim(instances)[1]
	belief <- infer.method(crf, ...)
	nll$value <- n.instances * belief$logZ - features$params %*% features$suff.stat
	g <- -features$suff.stat
	for (n in 1:crf$n.nodes)
	{
		for (s in 1:crf$n.states[n])
		{
			k <- features$node.pot[n,s,1]
			if (k) g[k] <- g[k] + n.instances * belief$node.bel[n,s]
		}
	}
	for (e in 1:crf$n.edges)
	{
		for (s1 in 1:crf$n.states[crf$edges[e,1]])
		{
			for (s2 in 1:crf$n.states[crf$edges[e,2]])
			{
				k <- features$edge.pot[[e]][s1,s2,1]
				if (k) g[k] <- g[k] + n.instances * belief$edge.bel[[e]][s1,s2]
			}
		}
	}
	nll$gradient <- g
	nll
}

get.potential <- function(crf, configuration)
	.Call("Get_Potential", crf, configuration)

get.logPotential <- function(crf, configuration)
	.Call("Get_LogPotential", crf, configuration)
