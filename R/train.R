make.features <- function(crf, n.nf = 1, n.ef = n.nf, n.par = 1)
{
	crf$n.nf <- n.nf
	crf$n.ef <- n.ef
	crf$node.par <- array(0, dim=c(crf$n.nodes, crf$max.state, n.nf))
	crf$edge.par <- lapply(1:crf$n.edges, function(i) array(0, dim=c(crf$n.states[crf$edges[i,1]], crf$n.states[crf$edges[i,2]], n.ef)))
	crf$n.par <- n.par
	crf$par <- numeric(crf$n.par)
	crf
}

make.node.pot <- function(crf, features = 1)
{
	features <- matrix(features, crf$n.nf, crf$n.nodes)
	index <- crf$node.par > 0
	node.par <- crf$node.par
	node.par[index] <- crf$par[node.par[index]]
	crf$node.pot <- exp(rowSums(node.par * as.vector(apply(features, 1, function(i) rep(i, crf$max.state))), dims=2))
	crf
}

make.edge.pot <- function(crf, features = 1)
{
	features <- matrix(features, crf$n.ef, crf$n.edges)
	make.edge.pot.i <- function(i)
	{
		index <- crf$edge.par[[i]] > 0
		edge.par <- crf$edge.par[[i]]
		edge.par[index] <- crf$par[edge.par[index]]
		edge.pot <- exp(rowSums(edge.par * rep(features[,i], each=length(crf$edge.pot[[i]])), dims=2))
	}
	crf$edge.pot <- lapply(1:crf$n.edges, make.edge.pot.i)
	crf
}

calc.frequency <- function(v, n)
	.Call("Calc_Frequency", v, n)

mrf.suff.stat <- function(crf, instances)
{
	suff.stat <- numeric(crf$n.par)
	n.instances <- dim(instances)[1]
	index <- cbind(rep(1:crf$n.nodes, each=n.instances), as.vector(instances), rep(1, length(instances)))
	index <- crf$node.par[index]
	suff.stat <- suff.stat + calc.frequency(index, crf$n.par)
	for (i in 1:crf$n.edges)
	{
		index <- cbind(instances[,crf$edges[i,1]], instances[,crf$edges[i,2]], rep(1, n.instances))
		index <- crf$edge.par[[i]][index]
		suff.stat <- suff.stat + calc.frequency(index, crf$n.par)
	}
	suff.stat
}

mrf.nll <- function(par, crf, instances, infer.method = infer.chain, ...)
{
	n.instances <- dim(instances)[1]
	crf$par <- par
	crf <- make.node.pot(crf)
	crf <- make.edge.pot(crf)
	belief <- infer.method(crf, ...)
	nll <- as.vector(n.instances * belief$logZ - par %*% crf$suff.stat)
	nll
}

mrf.gradient <- function(par, crf, instances, infer.method = infer.chain, ...)
{
	n.instances <- dim(instances)[1]
	crf$par <- par
	crf <- make.node.pot(crf)
	crf <- make.edge.pot(crf)
	belief <- infer.method(crf, ...)
	gradient <- -crf$suff.stat
	for (n in 1:crf$n.nodes)
	{
		for (s in 1:crf$n.states[n])
		{
			k <- crf$node.par[n,s,1]
			if (k) gradient[k] <- gradient[k] + n.instances * belief$node.bel[n,s]
		}
	}
	for (e in 1:crf$n.edges)
	{
		for (s1 in 1:crf$n.states[crf$edges[e,1]])
		{
			for (s2 in 1:crf$n.states[crf$edges[e,2]])
			{
				k <- crf$edge.par[[e]][s1,s2,1]
				if (k) gradient[k] <- gradient[k] + n.instances * belief$edge.bel[[e]][s1,s2]
			}
		}
	}
	gradient
}

crf.nll <- function(par, crf, instances, node.fea = 1, edge.fea = 1, infer.method = infer.chain, ...)
{
	n.instances <- dim(instances)[1]
	node.fea <- array(node.fea, dim=c(crf$n.nf, crf$n.nodes, n.instances))
	edge.fea <- array(edge.fea, dim=c(crf$n.ef, crf$n.edges, n.instances))
	crf$par <- par
	nll <- 0
	for (i in 1:n.instances)
	{
		crf <- make.node.pot(crf, node.fea[,,i])
		crf <- make.edge.pot(crf, edge.fea[,,i])
		belief <- infer.method(crf, ...)
		nll <- nll - get.logPotential(crf, instances[i,]) + belief$logZ
	}
	nll
}

crf.gradient <- function(par, crf, instances, node.fea = array(1, dim=c(crf$n.nf, crf$n.nodes, dim(instances)[1])), edge.fea = array(1, dim=c(crf$n.ef, crf$n.edges, dim(instances)[1])), infer.method = infer.chain, ...)
{
}

train.mrf <- function(crf, rain)
{
	solution <- optim(crf$par, mrf.nll, mrf.gradient, crf, rain, method = "L-BFGS-B")
	crf$par <- solution$par
	crf <- make.node.pot(crf)
	crf <- make.edge.pot(crf)
	crf
}

train.crf <- function(crf, rain)
{
}
