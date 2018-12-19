#include "CRF.h"

SEXP Infer_TRBP(SEXP _crf, SEXP _maxIter, SEXP _cutoff, SEXP _verbose, SEXP _maximize)
{
	int maxIter = INTEGER_POINTER(AS_INTEGER(_maxIter))[0];
	double cutoff = NUMERIC_POINTER(AS_NUMERIC(_cutoff))[0];
	int verbose = INTEGER_POINTER(AS_INTEGER(_verbose))[0];
	bool maximize = LOGICAL_POINTER(AS_LOGICAL(_maximize))[0];

	CRF crf(_crf);
	crf.Init_Belief();
	crf.Infer_TRBP(maxIter, cutoff, verbose, maximize);

	return(crf._belief);
}

void CRF::Infer_TRBP(int maxIter, double cutoff, int verbose, bool maximize)
{
	double *mu = (double *) R_alloc(nEdges, sizeof(double));
	double **scaleEdgePot = (double **) R_allocArray2<double>(nEdges, nEdgeStates);

	TRBP_Init(mu, scaleEdgePot);
	TRBP(mu, scaleEdgePot, maxIter, cutoff, verbose, maximize);
	TRBP_Messages2EdgeBel(mu, scaleEdgePot);
	TRBP_BetheFreeEnergy(mu);
}
