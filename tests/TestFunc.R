test.decode <- function(name, decode.method, crf, answer, ...)
{
	cat("  ", name, ": Decoding ... ", sep="")
	decode <- decode.method(crf, ...)

	if (all(decode == answer$decode)) {
		cat("Passed.\n")
	} else {
		cat("\n")
		stop("Decoding is incorrect!")
	}
}

test.infer <- function(name, infer.method, crf, answer, cutoff=1e-8, ...)
{
	cat("  ", name, ": Inferring ... ", sep="")
	belief <- infer.method(crf, ...)

	if (max(abs(c(belief$node.bel - answer$node.bel, belief$edge.bel - answer$edge.bel, belief$logZ - answer$logZ))) < cutoff) {
		cat("Passed.\n")
	} else {
		cat("\n")
		stop("Inference is incorrect!")
	}
}

test.sample <- function(name, sample.method, crf, answer, cutoff=0.01, size=10000, ...)
{
	cat("  ", name, ": Sampling ... ", sep="")
	samples <- sample.method(crf, size, ...)

	samples.node.bel <- array(0, dim=c(crf$n.nodes, crf$max.state))
	for (i in 1:crf$n.nodes)
		for (j in 1:crf$max.state)
			samples.node.bel[i,j] <- sum(samples[,i] == j)
	samples.node.bel = samples.node.bel / rowSums(samples.node.bel)
	if (mean(abs(samples.node.bel - answer$node.bel)) < cutoff) {
		cat("Passed.\n")
	} else {
		cat("\n")
		stop("Sampling may be incorrect!")
	}
}

save(test.decode, test.infer, test.sample, file="TestFunc.RData")
