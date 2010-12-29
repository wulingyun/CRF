infer.exact <- function(crf)
	.Call("Infer_Exact", crf)

infer.chain <- function(crf)
	.Call("Infer_Chain", crf)

infer.tree <- function(crf)
	.Call("Infer_Tree", crf)
