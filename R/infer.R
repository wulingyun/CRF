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
	newcrf <- clamp.crf(crf, clamped)
	b <- infer.method(newcrf, ...)
	belief$node.bel[newcrf$node.id, 1:newcrf$max.state] <- b$node.bel
	belief$edge.bel[1:newcrf$max.state, 1:newcrf$max.state, newcrf$edge.id] <- b$edge.bel
	belief$logZ <- b$logZ
	belief$node.bel[cbind(which(clamped != 0), clamped[clamped != 0])] <- 1
	e <- newcrf$original$edges
	e0 <- which(clamped[e[,1]] != 0 & clamped[e[,2]] != 0)
	e1 <- which(clamped[e[,1]] != 0 & clamped[e[,2]] == 0)
	e2 <- which(clamped[e[,1]] == 0 & clamped[e[,2]] != 0)
	if (length(e0) > 0) belief$edge.bel[cbind(clamped[e[e0,1]], clamped[e[e0,2]], e0)] <- 1
	if (length(e1) > 0) belief$edge.bel[cbind(rep(clamped[e[e1,1]], each=newcrf$max.state), 1:newcrf$max.state, rep(e1, each=newcrf$max.state))] <- t(b$node.bel[newcrf$node.map[e[e1,2]],])
	if (length(e2) > 0) belief$edge.bel[cbind(1:newcrf$max.state, rep(clamped[e[e2,2]], each=newcrf$max.state), rep(e2, each=newcrf$max.state))] <- t(b$node.bel[newcrf$node.map[e[e2,1]],])
	belief
}

infer.cutset <- function(crf, cutset, engine = "default")
{
	engine.id <- c("default"=-1, "none"=0, "exact"=1, "chain"=2, "tree"=3);
	clamped <- rep(0, crf$n.nodes)
	clamped[cutset] <- 1
	newcrf <- clamp.crf(crf, clamped)
	.Call("Infer_Cutset", newcrf, engine.id[engine])
}

infer.sample <- function(crf, sample.method, ...)
	.Call("Infer_Sample", crf, sample.method(crf, ...))

infer.lbp <- function(crf, max.iter = 10000, cutoff = 1e-4, verbose = 0)
	.Call("Infer_LBP", crf, max.iter, cutoff, verbose)

infer.trbp <- function(crf, max.iter = 10000, cutoff = 1e-4, verbose = 0)
	.Call("Infer_TRBP", crf, max.iter, cutoff, verbose)
