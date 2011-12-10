#include "CRF.h"

SEXP Decode_TRBP(SEXP _crf, SEXP _maxIter, SEXP _cutoff, SEXP _verbose)
{
	int maxIter = INTEGER_POINTER(AS_INTEGER(_maxIter))[0];
	double cutoff = NUMERIC_POINTER(AS_NUMERIC(_cutoff))[0];
	int verbose = INTEGER_POINTER(AS_INTEGER(_verbose))[0];

	CRF crf(_crf);
	crf.Init_Labels();
	crf.Init_NodeBel();
	crf.Decode_TRBP(maxIter, cutoff, verbose);

	return(crf._labels);
}

void CRF::Decode_TRBP(int maxIter, double cutoff, int verbose)
{
	double *mu = (double *) R_alloc(nEdges, sizeof(double));
	double **scaleEdgePot = (double **) allocArray2<double>(nEdges, nEdgeStates);

	TRBP_Weights(mu);
	TRBP_ScaleEdgePot(mu, scaleEdgePot);
	MessagesInit();
	TRBP(mu, scaleEdgePot, maxIter, cutoff, verbose, true);
	TRBP_Messages2NodeBel(mu);
	MaxOfMarginals();
}
