make.features <- function(crf, n.nf = 1, n.ef = 1)
{
	crf$n.nf <- n.nf
	crf$n.ef <- n.ef
	crf$node.par <- list()
	if (n.nf > 0) crf$node.par[[1]] <- matrix(0, crf$n.nodes, crf$max.state)
	if (n.nf > 1) crf$node.par[2:n.nf] <- crf$node.par[1]
	crf$edge.par <- list()
	if (n.ef > 0) crf$edge.par[[1]] <- lapply(1:crf$n.edges, function(i) matrix(0, crf$n.states[crf$edges[i,1]], crf$n.states[crf$edges[i,2]]))
	if (n.ef > 1) crf$edge.par[2:n.ef] <- crf$edge.par[1]
	crf
}

make.par <- function(crf, n.par = 1)
{
	crf$n.par <- n.par
	crf$par <- numeric(crf$n.par)
	crf$nll <- numeric(1)
	crf$gradient <- numeric(crf$n.par)
	crf
}

mrf.update <- function(crf)
	.Call("MRF_Update", crf)

crf.update <- function(crf, node.fea = NaN, edge.fea = NaN, node.ext = NaN, edge.ext = NaN)
	.Call("CRF_Update", crf, node.fea, edge.fea, node.ext, edge.ext)

mrf.stat <- function(crf, instances)
	.Call("MRF_Stat", crf, instances)

mrf.nll <- function(par, crf, instances, infer.method = infer.chain, ...)
	.Call("MRF_NLL", crf, par, instances, quote(infer.method(crf, ...)), environment())

crf.nll <- function(par, crf, instances, node.fea = NaN, edge.fea = NaN, node.ext = NaN, edge.ext = NaN, infer.method = infer.chain, ...)
	.Call("CRF_NLL", crf, par, instances, node.fea, edge.fea, node.ext, edge.ext, quote(infer.method(crf, ...)), environment())

gradient <- function(par, crf, ...)
	crf$gradient

train.mrf <- function(crf, instances, trace = 0)
{
	crf$par.stat <- mrf.stat(crf, instances)
	solution <- optim(crf$par, mrf.nll, gradient, crf, instances, method = "L-BFGS-B", control = list(trace = trace))
	crf$par <- solution$par
	mrf.update(crf)
	crf
}

train.crf <- function(crf, instances, node.fea = NaN, edge.fea = NaN, node.ext = NaN, edge.ext = NaN, trace = 0)
{
	solution <- optim(crf$par, crf.nll, gradient, crf, instances, node.fea, edge.fea, node.ext, edge.ext, method = "L-BFGS-B", control = list(trace = trace))
	crf$par <- solution$par
	crf.update(crf, node.fea[[1]], edge.fea[[1]], node.ext[[1]], edge.ext[[1]])
	crf
}
