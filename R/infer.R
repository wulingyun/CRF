infer.exact <- function(crf)
	.Call("Infer_Exact", crf)

infer.chain <- function(crf)
	.Call("Infer_Chain", crf)

infer.tree <- function(crf)
	.Call("Infer_Tree", crf)

infer.lbp <- function(crf, max.iter = 10000, cutoff = 1e-4, verbose = 0)
	.Call("Infer_LBP", crf, max.iter, cutoff, verbose)

infer.sample <- function(crf, sample.method, ...)
{
	s <- sample.method(crf, ...)
	belief <- list()
	belief$node.bel <- array(0, dim=c(crf$n.nodes, crf$max.state))
	for (i in 1:crf$n.nodes)
	{
		x <- table(s[,i])
		belief$node.bel[i, 1:length(x)] <- x / dim(s)[1]
	}
	belief$edge.bel <- array(0, dim=c(crf$max.state, crf$max.state, crf$n.edges))
	for (i in 1:crf$n.edges)
	{
		x <- table(s[,crf$edges[i,1]], s[,crf$edges[i,2]])
		belief$edge.bel[1:dim(x)[1], 1:dim(x)[2], i] <- x / dim(s)[1]
	}
	belief$logZ <- log(sum(apply(unique(s), 1, function(x) potential(crf, x))))
	belief
}
