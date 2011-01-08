infer.exact <- function(crf)
	.Call("Infer_Exact", crf)

infer.chain <- function(crf)
	.Call("Infer_Chain", crf)

infer.tree <- function(crf)
	.Call("Infer_Tree", crf)

infer.lbp <- function(crf, max.iter = 10000, cutoff = 1e-4, debug = 0)
	.Call("Infer_LBP", crf, max.iter, cutoff, debug)
