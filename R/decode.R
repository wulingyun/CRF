decode.exact <- function(crf)
	.Call("Decode_Exact", crf)

decode.chain <- function(crf)
	.Call("Decode_Chain", crf)

decode.tree <- function(crf)
	.Call("Decode_Tree", crf)

decode.lbp <- function(crf, max.iter = 10000, cutoff = 1e-4, verbose = 0)
	.Call("Decode_LBP", crf, max.iter, cutoff, verbose)

decode.infer <- function(crf, infer, ...)
	apply(infer(crf, ...)$node.bel, 1, which.max)

decode.sample <- function(crf, sample, ...)
{
	s <- sample(crf, ...)
	s[which.max(apply(s, 1, function(x) log.potential(crf, x))),]
}
