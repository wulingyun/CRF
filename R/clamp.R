clamp.crf <- function(crf, clamped)
{
	data <- new.env()
	if (!is.vector(clamped) || length(clamped) != crf$n.nodes)
		stop("'clamped' should be a vector of length ", crf$n.nodes, "!")
	if (any(clamped > crf$n.states | clamped < 0))
		stop("'clamped' has invalid value(s)!")

	data$original <- crf
	data$clamped <- clamped

	data$node.id <- which(clamped == 0)
	data$n.nodes <- length(data$node.id)
	data$node.map <- rep(0, crf$n.nodes)
	data$node.map[data$node.id] <- 1:data$n.nodes

	data$edge.id <- which(clamped[crf$edges[,1]] == 0 & clamped[crf$edges[,2]] == 0)
	data$n.edges <- length(data$edge.id)
	data$edges <- matrix(data$node.map[crf$edges[data$edge.id,]], ncol=2)
	data$edge.map <- rep(0, crf$n.edges)
	data$edge.map[data$edge.id] <- 1:data$n.edges

	.Call("Make_AdjInfo", data)

	data$n.states <- crf$n.states[data$node.id]
	data$max.state <- max(data$n.states)

	data$node.pot <- crf$node.pot[data$node.id, 1:data$max.state]
	data$edge.pot <- crf$edge.pot[data$edge.id]
	.Call("Clamp_Reset", data)

	class(data) <- c("CRF.clamped", "CRF")
	data
}

clamp.reset <- function(crf, clamped)
{
	if (is.na(class(crf)[1]) || class(crf)[1] != "CRF.clamped")
		stop("'crf' is not class CRF.clamped!")
	if (sum(xor(crf$clamped == 0, clamped == 0)) != 0)
		stop("'clamped' has different clamped structure!")
	if (any(clamped > crf$original$n.states | clamped < 0))
		stop("'clamped' has invalid clamped value(s)!")
	crf$clamped <- clamped
	.Call("Clamp_Reset", crf)
	crf
}

sub.crf <- function(crf, subset)
{
	data <- new.env()
	data$original <- crf

	data$node.id <- intersect(1:crf$n.nodes, unique(subset))
	data$n.nodes <- length(data$node.id)
	data$node.map <- rep(0, crf$n.nodes)
	data$node.map[data$node.id] <- 1:data$n.nodes

	data$edge.id <- which(data$node.map[crf$edges[,1]] != 0 & data$node.map[crf$edges[,2]] != 0)
	data$n.edges <- length(data$edge.id)
	data$edges <- matrix(data$node.map[crf$edges[data$edge.id,]], ncol=2)
	data$edge.map <- rep(0, crf$n.edges)
	data$edge.map[data$edge.id] <- 1:data$n.edges

	adj.info <- .Call("Make_AdjInfo", data)
	data$n.adj <- adj.info$n.adj
	data$adj.nodes <- adj.info$adj.nodes
	data$adj.edges <- adj.info$adj.edges

	data$n.states <- crf$n.states[data$node.id]
	data$max.state <- max(data$n.states)

	data$node.pot <- array(crf$node.pot[data$node.id, 1:data$max.state], dim=c(data$n.nodes, data$max.state))
	data$edge.pot <- crf$edge.pot[data$edge.id]

	class(data) <- c("CRF.sub", "CRF")
	data
}
