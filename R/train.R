make.features <- function(crf, n.nf = 1, n.ef = n.nf, n.par = 1)
{
	crf$n.nf <- n.nf
	crf$n.ef <- n.ef
	crf$node.par <- array(0, dim=c(crf$n.nodes, crf$max.state, n.nf))
	crf$edge.par <- lapply(1:crf$n.edges, function(i) array(0, dim=c(crf$n.states[crf$edges[i,1]], crf$n.states[crf$edges[i,2]], n.ef)))
	crf$n.par <- n.par
	crf$par <- numeric(crf$n.par)
	crf$w.par <- numeric(crf$n.par)
	crf
}

update.pot <- function(crf, node.fea = 1, edge.fea = 1)
{
	node.fea <- matrix(node.fea, crf$n.nf, crf$n.nodes)
	edge.fea <- matrix(edge.fea, crf$n.ef, crf$n.edges)
	.Call("Update_Pot", crf, node.fea, edge.fea)
}

update.par.stat <- function(crf, instances, node.fea = 1, edge.fea = 1)
{
	n.instances <- dim(instances)[1]
	node.fea <- array(node.fea, dim=c(crf$n.nf, crf$n.nodes, n.instances))
	edge.fea <- array(edge.fea, dim=c(crf$n.ef, crf$n.edges, n.instances))
	.Call("Update_ParStat", crf, n.instances, instances, node.fea, edge.fea)
}

mrf.nll <- function(par, crf, instances, infer.method = infer.chain, ...)
{
	n.instances <- dim(instances)[1]
	crf$par <- par
	update.pot(crf)
	belief <- infer.method(crf, ...)
	nll <- as.vector(n.instances * belief$logZ - par %*% crf$w.par)
	nll
}

mrf.gradient <- function(par, crf, instances, infer.method = infer.chain, ...)
{
	n.instances <- dim(instances)[1]
	crf$par <- par
	update.pot(crf)
	belief <- infer.method(crf, ...)
	gradient <- -crf$w.par
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
		update.pot(crf, node.fea[,,i], edge.fea[,,i])
		belief <- infer.method(crf, ...)
		nll <- nll - get.logPotential(crf, instances[i,]) + belief$logZ
	}
	nll
}

crf.gradient <- function(par, crf, instances, node.fea = 1, edge.fea = 1, infer.method = infer.chain, ...)
{
	n.instances <- dim(instances)[1]
	node.fea <- array(node.fea, dim=c(crf$n.nf, crf$n.nodes, n.instances))
	edge.fea <- array(edge.fea, dim=c(crf$n.ef, crf$n.edges, n.instances))
	crf$par <- par
	gradient <- -crf$w.par
	for (i in 1:n.instances)
	{
		update.pot(crf, node.fea[,,i], edge.fea[,,i])
		belief <- infer.method(crf, ...)
		for (n in 1:crf$n.nodes)
		{
			for (s in 1:crf$n.states[n])
			{
				k <- crf$node.par[n,s,1]
				if (k) gradient[k] <- gradient[k] + belief$node.bel[n,s]
			}
		}
		for (e in 1:crf$n.edges)
		{
			for (s1 in 1:crf$n.states[crf$edges[e,1]])
			{
				for (s2 in 1:crf$n.states[crf$edges[e,2]])
				{
					k <- crf$edge.par[[e]][s1,s2,1]
					if (k) gradient[k] <- gradient[k] + belief$edge.bel[[e]][s1,s2]
				}
			}
		}
	}
	gradient
}

train.mrf <- function(crf, instances)
{
	update.par.stat(crf, instances)
	solution <- optim(crf$par, mrf.nll, mrf.gradient, crf, instances, method = "L-BFGS-B")
	crf$par <- solution$par
	update.pot(crf)
	crf
}

train.crf <- function(crf, instances, node.fea = 1, edge.fea = 1)
{
	update.par.stat(crf, instances, node.fea, edge.fea)
	solution <- optim(crf$par, crf.nll, crf.gradient, crf, instances, node.fea, edge.fea, method = "L-BFGS-B")
	crf$par <- solution$par
	update.pot(crf)
	crf
}
