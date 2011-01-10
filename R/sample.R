sample.exact <- function(crf, size)
	.Call("Sample_Exact", crf, size)

sample.chain <- function(crf, size)
	.Call("Sample_Chain", crf, size)

sample.tree <- function(crf, size)
	.Call("Sample_Tree", crf, size)

sample.gibbs <- function(crf, size, burn.in = 1000, start = apply(crf$node.pot, 1, which.max))
	.Call("Sample_Gibbs", crf, size, burn.in, start)
