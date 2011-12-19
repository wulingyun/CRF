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
	belief$edge.bel <- lapply(1:crf$n.edges, function(i) array(0, dim=c(crf$n.states[crf$edges[i,1]], crf$n.states[crf$edges[i,2]])))
	newcrf <- clamp.crf(crf, clamped)
	b <- infer.method(newcrf, ...)
	belief$node.bel[newcrf$node.id, 1:newcrf$max.state] <- b$node.bel
	belief$edge.bel[newcrf$edge.id] <- b$edge.bel
	belief$logZ <- b$logZ
	belief$node.bel[cbind(which(clamped != 0), clamped[clamped != 0])] <- 1
	e <- newcrf$original$edges
	e0 <- which(clamped[e[,1]] != 0 & clamped[e[,2]] != 0)
	e1 <- which(clamped[e[,1]] != 0 & clamped[e[,2]] == 0)
	e2 <- which(clamped[e[,1]] == 0 & clamped[e[,2]] != 0)
	for (i in e0) belief$edge.bel[[i]][clamped[e[i,1]], clamped[e[i,2]]] <- 1
	for (i in e1) belief$edge.bel[[i]][clamped[e[i,1]],] <- b$node.bel[newcrf$node.map[e[i,2]], 1:crf$n.states[e[i,2]]]
	for (i in e2) belief$edge.bel[[i]][,clamped[e[i,2]]] <- b$node.bel[newcrf$node.map[e[i,1]], 1:crf$n.states[e[i,1]]]
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

infer.junction <- function(crf)
	.Call("Infer_Junction", crf)

infer.sample <- function(crf, sample.method, ...)
	.Call("Infer_Sample", crf, sample.method(crf, ...))

infer.lbp <- function(crf, max.iter = 10000, cutoff = 1e-4, verbose = 0)
	.Call("Infer_LBP", crf, max.iter, cutoff, verbose)

infer.trbp <- function(crf, max.iter = 10000, cutoff = 1e-4, verbose = 0)
	.Call("Infer_TRBP", crf, max.iter, cutoff, verbose)
