infer.exact <- function(crf)
	.Call("Infer_Exact", crf)

infer.chain <- function(crf)
	.Call("Infer_Chain", crf)

infer.tree <- function(crf)
	.Call("Infer_Tree", crf)

infer.conditional <- function(crf, clamped, infer.method, ...)
{
	belief <- list()
	belief$node.bel <- array(0, dim=c(crf$n.nodes, crf$max.state))
	belief$edge.bel <- array(0, dim=c(crf$max.state, crf$max.state, crf$n.edges))
	crf <- clamp.crf(crf, clamped)
	b <- infer.method(crf, ...)
	belief$node.bel[crf$node.id, 1:crf$max.state] <- b$node.bel
	belief$edge.bel[1:crf$max.state, 1:crf$max.state, crf$edge.id] <- b$edge.bel
	belief$logZ <- b$logZ
	belief$node.bel[cbind(which(clamped != 0), clamped[clamped != 0])] <- 1
	e <- crf$original$edges
	e0 <- which(clamped[e[,1]] != 0 & clamped[e[,2]] != 0)
	e1 <- which(clamped[e[,1]] != 0 & clamped[e[,2]] == 0)
	e2 <- which(clamped[e[,1]] == 0 & clamped[e[,2]] != 0)
	if (length(e0) > 0) belief$edge.bel[cbind(clamped[e[e0,1]], clamped[e[e0,2]], e0)] <- 1
	if (length(e1) > 0) belief$edge.bel[cbind(rep(clamped[e[e1,1]], each=crf$max.state), 1:crf$max.state, rep(e1, each=crf$max.state))] <- t(b$node.bel[crf$node.map[e[e1,2]],])
	if (length(e2) > 0) belief$edge.bel[cbind(1:crf$max.state, rep(clamped[e[e2,2]], each=crf$max.state), rep(e2, each=crf$max.state))] <- t(b$node.bel[crf$node.map[e[e2,1]],])
	belief
}

infer.cutset <- function(crf, cutset)
{
	clamped <- rep(0, crf$n.nodes)
	clamped[cutset] <- 1
	crf <- clamp.crf(crf, clamped)
	.Call("Infer_Cutset", crf)
}

infer.sample <- function(crf, sample.method, ...)
	.Call("Infer_Sample", crf, sample.method(crf, ...))

infer.lbp <- function(crf, max.iter = 10000, cutoff = 1e-4, verbose = 0)
	.Call("Infer_LBP", crf, max.iter, cutoff, verbose)

infer.trbp <- function(crf, max.iter = 10000, cutoff = 1e-4, verbose = 0)
	.Call("Infer_TRBP", crf, max.iter, cutoff, verbose)
