#' Make CRF features
#' 
#' Make the data structure of features
#' 
#' Make the data structure of features need for modeling and training 
#' 
#' @param crf
#' @param n.nf
#' @param n.ef
#' @return This function will return the same CRF.
#' 
#' @export
make.features <- function(crf, n.nf = 1, n.ef = 1)
{
	crf$n.nf <- n.nf
	crf$n.ef <- n.ef
	crf$node.par <- array(0, dim=c(crf$n.nodes, crf$max.state, n.nf))
	crf$edge.par <- lapply(1:crf$n.edges, function(i) array(0, dim=c(crf$n.states[crf$edges[i,1]], crf$n.states[crf$edges[i,2]], n.ef)))
	crf
}



#' Make CRF parameters
#' 
#' Make the data structure of parameters
#' 
#' Make the data structure of parameters need for modeling and training 
#' 
#' @param crf
#' @param n.par
#' @return This function will return the same CRF.
#' 
#' @export
make.par <- function(crf, n.par = 1)
{
	crf$n.par <- n.par
	crf$par <- numeric(crf$n.par)
	crf$nll <- numeric(1)
	crf$gradient <- numeric(crf$n.par)
	crf
}



#' Update MRF potentials
#' 
#' Update node.pot and edge.pot of MRF model 
#' 
#' Update node.pot and edge.pot of MRF model 
#' 
#' @param crf
#' @return This function will directly modify the CRF. Do not use the returned value.
#' 
#' @export
mrf.update <- function(crf)
	.Call(MRF_Update, crf)



#' Update CRF potentials
#' 
#' Update node.pot and edge.pot of CRF model 
#' 
#' Update node.pot and edge.pot of CRF model 
#' 
#' @param crf
#' @param node.fea
#' @param edge.fea
#' @param node.ext
#' @param edge.ext
#' @return This function will directly modify the CRF. Do not use the returned value.
#' 
#' @export
crf.update <- function(crf, node.fea = NaN, edge.fea = NaN, node.ext = NaN, edge.ext = NaN)
	.Call(CRF_Update, crf, node.fea, edge.fea, node.ext, edge.ext)



#' Calculate MRF sufficient statistics
#' 
#' Calculate the sufficient statistics of MRF model 
#' 
#' Calculate the sufficient statistics of MRF model 
#' 
#' @param crf
#' @param instances
#' @return This function will return the value of MRF sufficient statistics.
#' 
#' @export
mrf.stat <- function(crf, instances)
	.Call(MRF_Stat, crf, instances)



#' Calculate MRF negative log-likelihood
#' 
#' Calculate the negative log-likelihood of MRF model 
#' 
#' Calculate the negative log-likelihood of MRF model 
#' 
#' @param crf
#' @param par
#' @param instances
#' @param infer.method
#' @param ...
#' @return This function will return the value of MRF negative log-likilihood.
#' 
#' @export
mrf.nll <- function(par, crf, instances, infer.method = infer.chain, ...)
	.Call(MRF_NLL, crf, par, instances, quote(infer.method(crf, ...)), environment())



#' Calculate CRF negative log likelihood
#' 
#' Calculate the negative log likelihood of CRF model
#' 
#' Calculate the negative log likelihood of CRF model
#' 
#' @param crf
#' @param par
#' @param instances
#' @param node.fea
#' @param edge.fea
#' @param node.ext
#' @param edge.ext
#' @param infer.method
#' @param ...
#' @return This function will return the value of CRF negative log-likelihood.
#' 
#' @export
crf.nll <- function(par, crf, instances, node.fea = NaN, edge.fea = NaN, node.ext = NaN, edge.ext = NaN, infer.method = infer.chain, ...)
	.Call(CRF_NLL, crf, par, instances, node.fea, edge.fea, node.ext, edge.ext, quote(infer.method(crf, ...)), environment())



#' Calculate CRF negative log-likelihood gradient
#' 
#' Calculate the gradient of negative log likelihood of CRF model
#' 
#' Calculate the gradient of negative log likelihood of CRF model. 
#' This function is used by optimization algorithm in training.
#' 
#' @param par
#' @param crf
#' @param ...
#' @return This function will return the gradient of CRF negative log-likelihood.
#' 
#' @export
gradient <- function(par, crf, ...)
	crf$gradient



#' Train MRF model
#' 
#' Train the MRF model to estimate the parameters
#' 
#' This function trains the Markov Random Fields (MRF) model, which is a simple variant of CRF model. 
#' 
#' @param crf
#' @param instances
#' @param trace
#' @return This function will return the same CRF.
#' 
#' @export
train.mrf <- function(crf, instances, trace = 0)
{
	crf$par.stat <- mrf.stat(crf, instances)
	solution <- optim(crf$par, mrf.nll, gradient, crf, instances, method = "L-BFGS-B", control = list(trace = trace))
	crf$par <- solution$par
	mrf.update(crf)
	crf
}



#' Train CRF model
#' 
#' Train the CRF model to estimate the parameters
#' 
#' This function train the CRF model. 
#' 
#' @param crf
#' @param instances
#' @param trace
#' @param node.fea
#' @param edge.fea
#' @param node.ext
#' @param edge.ext
#' @return This function will return the same CRF.
#' 
#' @export
train.crf <- function(crf, instances, node.fea = NaN, edge.fea = NaN, node.ext = NaN, edge.ext = NaN, trace = 0)
{
	solution <- optim(crf$par, crf.nll, gradient, crf, instances, node.fea, edge.fea, node.ext, edge.ext, method = "L-BFGS-B", control = list(trace = trace))
	crf$par <- solution$par
	crf.update(crf, node.fea[[1]], edge.fea[[1]], node.ext[[1]], edge.ext[[1]])
	crf
}
