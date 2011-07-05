sample.exact <- function(crf, size)
	.Call("Sample_Exact", crf, size)

sample.chain <- function(crf, size)
	.Call("Sample_Chain", crf, size)

sample.tree <- function(crf, size)
	.Call("Sample_Tree", crf, size)

sample.conditional <- function(crf, size, clamped, sample.method, ...)
{
	newcrf <- clamp.crf(crf, clamped)
	s <- sample.method(newcrf, size, ...)
	samples <- matrix(rep(clamped, nrow(s)), nrow=nrow(s), ncol=length(clamped), byrow=TRUE)
	samples[,newcrf$node.id] <- s
	samples
}

sample.cutset <- function(crf, size, cutset, engine = "default")
{
	engine.id <- c("default"=-1, "none"=0, "exact"=1, "chain"=2, "tree"=3);
	clamped <- rep(0, crf$n.nodes)
	clamped[cutset] <- 1
	newcrf <- clamp.crf(crf, clamped)
	.Call("Sample_Cutset", newcrf, size, engine.id[engine])
}

sample.gibbs <- function(crf, size, burn.in = 1000, start = apply(crf$node.pot, 1, which.max))
	.Call("Sample_Gibbs", crf, size, burn.in, start)
