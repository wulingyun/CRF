library(CRF)

nNodes <- 4
nStates <- 2

adj <- matrix(nrow=nNodes, ncol=nNodes)
adj[1,2] <- 1
adj[2,1] <- 1
adj[2,3] <- 1
adj[3,2] <- 1
adj[3,4] <- 1
adj[4,3] <- 1

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

decode <- decode.exact(crf)
print("Optimal decoding")
print(decode)

belief <- infer.exact(crf)
print(belief)

samples <- sample.exact(crf, 10000)
print(samples)

print(c(sum(samples[,1] == 1), sum(samples[,1] == 2)))
print(c(sum(samples[,2] == 1), sum(samples[,2] == 2)))
print(c(sum(samples[,3] == 1), sum(samples[,3] == 2)))
print(c(sum(samples[,4] == 1), sum(samples[,4] == 2)))
