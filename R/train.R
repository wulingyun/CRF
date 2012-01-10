make.features <- function(crf, n.nf = 1, n.ef = n.nf, n.par = 1)
{
	crf$n.nf <- n.nf
	crf$n.ef <- n.ef
	crf$node.fea <- array(0, dim=c(crf$n.nodes, crf$max.state, n.nf))
	crf$edge.fea <- lapply(1:crf$n.edges, function(i) array(0, dim=c(crf$n.states[crf$edges[i,1]], crf$n.states[crf$edges[i,2]], n.ef)))
	crf$n.par <- n.par
	crf$par <- numeric(crf$n.par)
	crf
}

make.node.pot <- function(crf, instance = matrix(1, crf$n.nodes, crf$n.nf))
{
	index <- crf$node.fea > 0
	feature <- crf$node.fea
	feature[index] <- crf$par[feature[index]]
	crf$node.pot <- exp(rowSums(feature * as.vector(apply(instance, 2, function(i) rep(i, crf$max.state))), dims=2))
	crf
}

make.edge.pot <- function(crf, instance = matrix(1, crf$n.edges, crf$n.ef))
{
	make.edge.pot.i <- function(i)
	{
		index <- crf$edge.fea[[i]] > 0
		feature <- crf$edge.fea[[i]]
		feature[index] <- crf$par[feature[index]]
		edge.pot <- exp(rowSums(feature * rep(instance[i,], each=length(crf$edge.pot[[i]])), dims=2))
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
	index <- crf$node.fea[index]
	suff.stat <- suff.stat + calc.frequency(index, crf$n.par)
	for (i in 1:crf$n.edges)
	{
		index <- cbind(instances[,crf$edges[i,1]], instances[,crf$edges[i,2]], rep(1, n.instances))
		index <- crf$edge.fea[[i]][index]
		suff.stat <- suff.stat + calc.frequency(index, crf$n.par)
	}
	suff.stat
}

mrf.nll <- function(par, crf, instances, infer.method=infer.chain, ...)
{
	n.instances <- dim(instances)[1]
	crf$par <- par
	crf <- make.node.pot(crf)
	crf <- make.edge.pot(crf)
	belief <- infer.method(crf, ...)
	nll <- as.vector(n.instances * belief$logZ - par %*% crf$suff.stat)
	nll
}

mrf.nll.gradient <- function(par, crf, instances, infer.method=infer.chain, ...)
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
			k <- crf$node.fea[n,s,1]
			if (k) gradient[k] <- gradient[k] + n.instances * belief$node.bel[n,s]
		}
	}
	for (e in 1:crf$n.edges)
	{
		for (s1 in 1:crf$n.states[crf$edges[e,1]])
		{
			for (s2 in 1:crf$n.states[crf$edges[e,2]])
			{
				k <- crf$edge.fea[[e]][s1,s2,1]
				if (k) gradient[k] <- gradient[k] + n.instances * belief$edge.bel[[e]][s1,s2]
			}
		}
	}
	gradient
}

train.mrf <- function(crf, rain)
{
	solution <- optim(crf$par, mrf.nll, mrf.nll.gradient, crf, rain, method = "L-BFGS-B")
	crf$par <- solution$par
	crf <- make.node.pot(crf)
	crf <- make.edge.pot(crf)
	crf
}
