decode.exact <- function(crf)
	.Call("Decode_Exact", crf)

decode.chain <- function(crf)
	.Call("Decode_Chain", crf)

decode.tree <- function(crf)
	.Call("Decode_Tree", crf)

decode.conditional <- function(crf, clamped, decode.method, ...)
{
	newcrf <- clamp.crf(crf, clamped)
	decode <- clamped
	decode[newcrf$node.id] <- decode.method(newcrf, ...)
	decode
}

decode.cutset <- function(crf, cutset, engine = "default", start = apply(crf$node.pot, 1, which.max))
{
	engine.id <- c("default"=-1, "none"=0, "exact"=1, "chain"=2, "tree"=3);
	clamped <- rep(0, crf$n.nodes)
	clamped[cutset] <- 1
	newcrf <- clamp.crf(crf, clamped)
	.Call("Decode_Cutset", newcrf, engine.id[engine], start)
}

decode.junction <- function(crf)
	.Call("Decode_Junction", crf);

decode.sample <- function(crf, sample.method, ...)
	.Call("Decode_Sample", crf, sample.method(crf, ...))

decode.marginal <- function(crf, infer.method, ...)
	apply(infer.method(crf, ...)$node.bel, 1, which.max)

decode.lbp <- function(crf, max.iter = 10000, cutoff = 1e-4, verbose = 0)
	.Call("Decode_LBP", crf, max.iter, cutoff, verbose)

decode.trbp <- function(crf, max.iter = 10000, cutoff = 1e-4, verbose = 0)
	.Call("Decode_TRBP", crf, max.iter, cutoff, verbose)

decode.greedy <- function(crf, restart = 0, start = apply(crf$node.pot, 1, which.max))
	.Call("Decode_Greedy", crf, restart, start)

decode.icm <- function(crf, restart = 0, start = apply(crf$node.pot, 1, which.max))
	.Call("Decode_ICM", crf, restart, start)

decode.block <- function(crf, blocks, decode.method = decode.tree, restart = 0, start = apply(crf$node.pot, 1, which.max), ...)
{
	y <- integer(crf$n.nodes)
	y[] <- start

	newcrf <- list()
	for (i in 1:length(blocks))
	{
		blocks[[i]] <- sort(blocks[[i]])
		clamped <- y
		clamped[blocks[[i]]] <- 0
		newcrf[[i]] <- clamp.crf(crf, clamped)
	}

	maxPot <- -1
	decode <- y
	restart <- max(0, restart)
	for (iter in 0:restart)
	{
		done = F
		while(!done)
		{
			done = T
			for (i in 1:length(blocks))
			{
				clamped <- y
				clamped[blocks[[i]]] <- 0
				newcrf[[i]] <- clamp.reset(newcrf[[i]], clamped)
				temp <- decode.method(newcrf[[i]], ...)
				if (any(temp != y[blocks[[i]]]))
				{
					y[blocks[[i]]] <- temp
					done = F
				}
			}
		}

		pot <- get.potential(crf, y)
		if (pot > maxPot)
		{
			maxPot <- pot
			decode <- y
		}

		if (iter < restart)
		{
			y <- ceiling(runif(crf$n.nodes) * crf$n.states)
		}
	}
	decode
}

decode.ilp <- function(crf, lp.rounding = FALSE)
{
	vmap.nodes <- matrix(nrow=crf$n.nodes, ncol=2)
	vmap.edges <- matrix(nrow=crf$n.edges, ncol=2)
	n <- 0
	for (i in 1:crf$n.nodes)
	{
		vmap.nodes[i, 1] <- n + 1
		n <- n + crf$n.states[i]
		vmap.nodes[i, 2] <- n
	}
	vnum.nodes <- n
	for (i in 1:crf$n.edges)
	{
		vmap.edges[i, 1] <- n + 1
		n <- n + crf$n.states[crf$edges[i,1]] * crf$n.states[crf$edges[i,2]]
		vmap.edges[i, 2] <- n
	}
	vnum.edges <- n - vnum.nodes
	vnum.total <- n

	obj <- rep(0, vnum.total)
	for (i in 1:crf$n.nodes)
	{
		obj[vmap.nodes[i,1]:vmap.nodes[i,2]] <- -log(crf$node.pot[i,1:crf$n.states[i]])
	}
	for (i in 1:crf$n.edges)
	{
		obj[vmap.edges[i,1]:vmap.edges[i,2]] <- -log(crf$edge.pot[[i]])
	}
	obj[is.infinite(obj)] <- 1000

	cnum.nodes <- crf$n.nodes
	cnum.edges <- sum(crf$n.states[crf$edges])
	cnum.total <- cnum.nodes + cnum.edges
	mat <- matrix(0, nrow=cnum.total, ncol=vnum.total)
	for (i in 1:crf$n.nodes)
	{
		mat[i, vmap.nodes[i,1]:vmap.nodes[i,2]] <- 1
	}
	n <- cnum.nodes
	for (i in 1:crf$n.edges)
	{
		n1 <- crf$edges[i,1]
		n2 <- crf$edges[i,2]
		for (j in 1:crf$n.states[n1])
		{
			n <- n + 1
			mat[n, vmap.nodes[n1,1]-1+j] <- -1
			mat[n, seq.int(vmap.edges[i,1]-1+j, vmap.edges[i,2], by=crf$n.states[n1])] <- 1
		}
		for (j in 1:crf$n.states[n2])
		{
			n <- n + 1
			mat[n, vmap.nodes[n2,1]-1+j] <- -1
			mat[n, (vmap.edges[i,1]+(j-1)*crf$n.states[n1]):(vmap.edges[i,1]-1+j*crf$n.states[n1])] <- 1
		}
	}

	dir <- rep("==", cnum.total)
	rhs <- rep(0, cnum.total)
	rhs[1:cnum.nodes] <- 1

	if (lp.rounding)
	{
		types <- rep("C", vnum.total)
		bounds <- list(upper = list(ind = 1:vnum.total, val = rep(1, vnum.total)))
	}
	else
	{
		types <- rep("B", vnum.total)
		bounds <- NULL
	}

	result <- Rglpk::Rglpk_solve_LP(obj, mat, dir, rhs, types = types, bounds = bounds)

	if (result$status != 0)
	{
		warning("LP solution is not optimal.")
	}

	node.bel <- matrix(0, nrow=crf$n.nodes, ncol=crf$max.state)
	for (i in 1:crf$n.nodes)
	{
		node.bel[i, 1:crf$n.states[i]] <- result$solution[vmap.nodes[i,1]:vmap.nodes[i,2]]
	}
	apply(node.bel, 1, which.max)
}
