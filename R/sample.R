sample.exact <- function(crf, size)
	.Call("Sample_Exact", crf, size)

sample.chain <- function(crf, size)
	.Call("Sample_Chain", crf, size)

sample.tree <- function(crf, size)
	.Call("Sample_Tree", crf, size)

sample.conditional <- function(crf, size, clamped, sample.method, ...)
{
	crf <- clamp.crf(crf, clamped)
	s <- sample.method(crf, size, ...)
	samples <- matrix(rep(clamped, nrow(s)), nrow=nrow(s), ncol=length(clamped), byrow=TRUE)
	samples[,crf$node.id] <- s
	samples
}

sample.cutset <- function(crf, size, cutset, is.chain = FALSE)
{
	clamped <- rep(0, crf$n.nodes)
	clamped[cutset] <- 1
	crf <- clamp.crf(crf, clamped)
	.Call("Sample_Cutset", crf, size, is.chain)
}

sample.cutsetChain <- function(crf, size, cutset)
{
	clamped <- rep(0, crf$n.nodes)
	clamped[cutset] <- 1
	crf <- clamp.crf(crf, clamped)
	.Call("Sample_CutsetChain", crf, size)
}

sample.gibbs <- function(crf, size, burn.in = 1000, start = apply(crf$node.pot, 1, which.max))
	.Call("Sample_Gibbs", crf, size, burn.in, start)
