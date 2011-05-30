library(CRF)
data(Chain)
crf <- chain.crf
answer <- chain.answer

print("Decoding ...")
decode <- decode.chain(crf)

if (all(decode == answer$decode)) {
	print("  Pass.")
} else {
	stop("Decoding is incorrect!")
}


print("Inferring ...")
belief <- infer.chain(crf)

if (max(abs(c(belief$node.bel - answer$node.bel, belief$edge.bel - answer$edge.bel, belief$logZ - answer$logZ))) < 1e-8) {
	print("  Pass.")
} else {
	stop("Inference is incorrect!")
}

print("Sampling ...")
samples <- sample.chain(crf, 10000)

samples.node.bel <- array(0, dim=c(crf$n.nodes, crf$max.state))
for (i in 1:crf$n.nodes)
	for (j in 1:crf$max.state)
		samples.node.bel[i,j] <- sum(samples[,i] == j)
samples.node.bel = samples.node.bel / rowSums(samples.node.bel)
if (mean(abs(samples.node.bel - answer$node.bel)) < 0.01) {
	print("  Pass.")
} else {
	stop("Sampling may be incorrect!")
}
