decode.exact <- function(crf)
	.Call("Decode_Exact", crf)

decode.chain <- function(crf)
	.Call("Decode_Chain", crf)

decode.tree <- function(crf)
	.Call("Decode_Tree", crf)

decode.lbp <- function(crf, max.iter = 10000, cutoff = 1e-4, verbose = 0)
	.Call("Decode_LBP", crf, max.iter, cutoff, verbose)

decode.marginal <- function(crf, infer.method, ...)
	apply(infer.method(crf, ...)$node.bel, 1, which.max)

decode.sample <- function(crf, sample.method, ...)
{
	s <- sample.method(crf, ...)
	s[which.max(apply(s, 1, function(x) log.potential(crf, x))),]
}
