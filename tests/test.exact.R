library(CRF)

nNodes <- 4
nStates <- 2

adj <- matrix(0, nrow=nNodes, ncol=nNodes)
for (i in 1:(nNodes-1))
{
	adj[i,i+1] <- 1
	adj[i+1,i] <- 1
}

crf <- make.crf(adj, nStates)

crf$node.pot[1,] <- c(1, 3)
crf$node.pot[2,] <- c(9, 1)
crf$node.pot[3,] <- c(1, 3)
crf$node.pot[4,] <- c(9, 1)

for (i in 1:crf$n.edges)
{
   crf$edge.pot[1,,i] <- c(2, 1)
   crf$edge.pot[2,,i] <- c(1, 2)
}

answer <-
structure(list(decode = c(2L, 1L, 1L, 1L), node.bel = structure(c(0.359630606860158, 
0.843007915567282, 0.486279683377309, 0.881002638522427, 0.640369393139842, 
0.156992084432718, 0.513720316622691, 0.118997361477573), .Dim = c(4L, 
2L)), edge.bel = structure(c(0.337203166226913, 0.505804749340369, 
0.0224274406332454, 0.134564643799472, 0.451187335092348, 0.0350923482849604, 
0.391820580474934, 0.121899736147757, 0.460686015831135, 0.420316622691293, 
0.0255936675461741, 0.0934036939313984), .Dim = c(2L, 2L, 3L)), 
    logZ = 8.24012129807647), .Names = c("decode", "node.bel", 
"edge.bel", "logZ"))

print("Decoding ...")
decode <- decode.exact(crf)

if (all(decode == answer$decode)) {
	print("  Pass.")
} else {
	stop("Decoding is incorrect!")
}

print("Inferring ...")
belief <- infer.exact(crf)

if (max(abs(c(belief$node.bel - answer$node.bel, belief$edge.bel - answer$edge.bel, belief$logZ - answer$logZ))) < 1e-8) {
	print("  Pass.")
} else {
	stop("Inference is incorrect!")
}

print("Sampling ...")
samples <- sample.exact(crf, 10000)

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
