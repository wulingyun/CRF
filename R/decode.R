decode.exact <- function(crf)
	.Call("Decode_Exact", crf)

decode.chain <- function(crf)
	.Call("Decode_Chain", crf)

decode.tree <- function(crf)
	.Call("Decode_Tree", crf)

decode.conditional <- function(crf, clamped, decode.method, ...)
{
	newcrf <- clamp.crf(crf, clamped)
	decode <- clamped
	decode[newcrf$node.id] <- decode.method(newcrf, ...)
	decode
}

decode.cutset <- function(crf, cutset, engine = "default", start = apply(crf$node.pot, 1, which.max))
{
	engine.id <- c("default"=-1, "none"=0, "exact"=1, "chain"=2, "tree"=3);
	clamped <- rep(0, crf$n.nodes)
	clamped[cutset] <- 1
	newcrf <- clamp.crf(crf, clamped)
	.Call("Decode_Cutset", newcrf, engine.id[engine], start)
}

decode.sample <- function(crf, sample.method, ...)
	.Call("Decode_Sample", crf, sample.method(crf, ...))

decode.marginal <- function(crf, infer.method, ...)
	apply(infer.method(crf, ...)$node.bel, 1, which.max)

decode.lbp <- function(crf, max.iter = 10000, cutoff = 1e-4, verbose = 0)
	.Call("Decode_LBP", crf, max.iter, cutoff, verbose)

decode.trbp <- function(crf, max.iter = 10000, cutoff = 1e-4, verbose = 0)
	.Call("Decode_TRBP", crf, max.iter, cutoff, verbose)

decode.greedy <- function(crf, restart = 0, start = apply(crf$node.pot, 1, which.max))
	.Call("Decode_Greedy", crf, restart, start)

decode.icm <- function(crf, restart = 0, start = apply(crf$node.pot, 1, which.max))
	.Call("Decode_ICM", crf, restart, start)
