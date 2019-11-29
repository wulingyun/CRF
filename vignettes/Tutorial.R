## ------------------------------------------------------------------------
library(CRF)

## ------------------------------------------------------------------------
n.nodes <- 10
n.states <- 2
prior.prob <- c(0.8, 0.2)
trans.prob <- matrix(0, nrow=2, ncol=2)
trans.prob[1,] <- c(0.95, 0.05)
trans.prob[2,] <- c(0.05, 0.95)

## ------------------------------------------------------------------------
prior.prob

## ------------------------------------------------------------------------
trans.prob

## ------------------------------------------------------------------------
adj <- matrix(0, n.nodes, n.nodes)
for (i in 1:(n.nodes-1))
{
	adj[i, i+1] <- 1
}

## ------------------------------------------------------------------------
adj

## ------------------------------------------------------------------------
mc <- make.crf(adj, n.states)

## ------------------------------------------------------------------------
mc$node.pot[1,] <- prior.prob
for (i in 1:mc$n.edges)
{
  mc$edge.pot[[i]] <- trans.prob
}

## ------------------------------------------------------------------------
mc.samples <- sample.chain(mc, 10000)
mc.samples[1:10, ]

## ------------------------------------------------------------------------
mrf.new <- make.crf(adj, n.states)

## ------------------------------------------------------------------------
mrf.new <- make.features(mrf.new)
mrf.new <- make.par(mrf.new, 4)

## ------------------------------------------------------------------------
mrf.new$node.par[1,1,1] <- 1
for (i in 1:mrf.new$n.edges)
{
	mrf.new$edge.par[[i]][1,1,1] <- 2
	mrf.new$edge.par[[i]][1,2,1] <- 3
	mrf.new$edge.par[[i]][2,1,1] <- 4
}

## ------------------------------------------------------------------------
mrf.new <- train.mrf(mrf.new, mc.samples)

## ------------------------------------------------------------------------
mrf.new$par

## ------------------------------------------------------------------------
mrf.new$node.pot <- mrf.new$node.pot / rowSums(mrf.new$node.pot)
mrf.new$edge.pot[[1]] <- mrf.new$edge.pot[[1]] / rowSums(mrf.new$edge.pot[[1]])

## ------------------------------------------------------------------------
mrf.new$node.pot[1,]

## ------------------------------------------------------------------------
mrf.new$edge.pot[[1]]

## ------------------------------------------------------------------------
emmis.prob <- matrix(0, nrow=2, ncol=4)
emmis.prob[1,] <- c(0.59, 0.25, 0.15, 0.01)
emmis.prob[2,] <- c(0.01, 0.15, 0.25, 0.59)
emmis.prob

## ------------------------------------------------------------------------
hmm.samples <- mc.samples
hmm.samples[mc.samples == 1] <- sample.int(4, sum(mc.samples == 1), replace = TRUE, prob=emmis.prob[1,])
hmm.samples[mc.samples == 2] <- sample.int(4, sum(mc.samples == 2), replace = TRUE, prob=emmis.prob[2,])
hmm.samples[1:10,]

## ------------------------------------------------------------------------
crf.new <- make.crf(adj, n.states)

## ------------------------------------------------------------------------
crf.new <- make.features(crf.new, 5, 1)
crf.new <- make.par(crf.new, 8)

## ------------------------------------------------------------------------
crf.new$node.par[1,1,1] <- 1
for (i in 1:crf.new$n.edges)
{
	crf.new$edge.par[[i]][1,1,] <- 2
	crf.new$edge.par[[i]][1,2,] <- 3
	crf.new$edge.par[[i]][2,1,] <- 4
}
crf.new$node.par[,1,2] <- 5
crf.new$node.par[,1,3] <- 6
crf.new$node.par[,1,4] <- 7
crf.new$node.par[,1,5] <- 8

## ------------------------------------------------------------------------
hmm.nf <- lapply(1:dim(hmm.samples)[1], function(i) matrix(1, crf.new$n.nf, crf.new$n.nodes))
for (i in 1:dim(hmm.samples)[1])
{
	hmm.nf[[i]][2, hmm.samples[i,] != 1] <- 0
	hmm.nf[[i]][3, hmm.samples[i,] != 2] <- 0
	hmm.nf[[i]][4, hmm.samples[i,] != 3] <- 0
	hmm.nf[[i]][5, hmm.samples[i,] != 4] <- 0
}
hmm.ef <- lapply(1:dim(hmm.samples)[1], function(i) matrix(1, crf.new$n.ef, crf.new$n.edges))

## ------------------------------------------------------------------------
crf.new <- train.crf(crf.new, mc.samples, hmm.nf, hmm.ef)

## ------------------------------------------------------------------------
crf.new$par

## ------------------------------------------------------------------------
hmm.infer <- matrix(0, nrow=dim(hmm.samples)[1], ncol=dim(hmm.samples)[2])
for (i in 1:dim(hmm.samples)[1])
{
  crf.new <- crf.update(crf.new, hmm.nf[[i]], hmm.ef[[i]])
  hmm.infer[i,] <- decode.chain(crf.new)
}

## ------------------------------------------------------------------------
sum(hmm.infer != mc.samples)

## ------------------------------------------------------------------------
crf.new <- train.crf(crf.new, mc.samples, hmm.nf, hmm.ef, infer.method = infer.lbp)

