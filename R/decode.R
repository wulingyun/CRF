decode.exact <- function(crf)
	.Call("Decode_Exact", crf)

decode.chain <- function(crf)
	.Call("Decode_Chain", crf)

decode.tree <- function(crf)
	.Call("Decode_Tree", crf)

decode.lbp <- function(crf, max.iter = 10000, cutoff = 1e-4, debug = 0)
	.Call("Decode_LBP", crf, max.iter, cutoff, debug)
