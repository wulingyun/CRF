sample.exact <- function(crf, size)
	.Call("Sample_Exact", crf, size)

sample.chain <- function(crf, size)
	.Call("Sample_Chain", crf, size)

sample.tree <- function(crf, size)
	.Call("Sample_Tree", crf, size)
