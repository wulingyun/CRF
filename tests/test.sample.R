library(CRF)
data(Small)
crf <- small.crf
answer <- small.answer

print("Decoding ...")
decode <- decode.sample(crf, sample.exact, 10000)

if (all(decode == answer$decode)) {
	print("  Pass.")
} else {
	stop("Decoding is incorrect!")
}

print("Inferring ...")
belief <- infer.sample(crf, sample.exact, 10000000)

if (max(abs(c(belief$node.bel - answer$node.bel, belief$edge.bel - answer$edge.bel, belief$logZ - answer$logZ))) < 1e-3) {
	print("  Pass.")
} else {
	stop("Inference is incorrect!")
}
