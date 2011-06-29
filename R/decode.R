decode.exact <- function(crf)
	.Call("Decode_Exact", crf)

decode.chain <- function(crf)
	.Call("Decode_Chain", crf)

decode.tree <- function(crf)
	.Call("Decode_Tree", crf)

decode.conditional <- function(crf, clamped, decode.method, ...)
{
	crf <- clamp.crf(crf, clamped)
	decode <- clamped
	decode[crf$node.id] <- decode.method(crf, ...)
	decode
}

decode.cutset <- function(crf, cutset)
{
	clamped <- rep(0, crf$n.nodes)
	clamped[cutset] <- 1
	crf <- clamp.crf(crf, clamped)
	.Call("Decode_Cutset", crf)
}

decode.cutsetChain <- function(crf, cutset)
{
	clamped <- rep(0, crf$n.nodes)
	clamped[cutset] <- 1
	crf <- clamp.crf(crf, clamped)
	.Call("Decode_CutsetChain", crf)
}

decode.sample <- function(crf, sample.method, ...)
	.Call("Decode_Sample", crf, sample.method(crf, ...))

decode.marginal <- function(crf, infer.method, ...)
	apply(infer.method(crf, ...)$node.bel, 1, which.max)

decode.lbp <- function(crf, max.iter = 10000, cutoff = 1e-4, verbose = 0)
	.Call("Decode_LBP", crf, max.iter, cutoff, verbose)

decode.trbp <- function(crf, max.iter = 10000, cutoff = 1e-4, verbose = 0)
	.Call("Decode_TRBP", crf, max.iter, cutoff, verbose)
