#include "CRF.h"

SEXP Decode_LBP(SEXP _crf, SEXP _maxIter, SEXP _cutoff, SEXP _verbose)
{
	int maxIter = INTEGER_POINTER(AS_INTEGER(_maxIter))[0];
	double cutoff = NUMERIC_POINTER(AS_NUMERIC(_cutoff))[0];
	int verbose = INTEGER_POINTER(AS_INTEGER(_verbose))[0];

	CRF crf(_crf);
	crf.Init_Labels();
	crf.Init_NodeBel();
	crf.Decode_LBP(maxIter, cutoff, verbose);

	return(crf._labels);
}

void CRF::Decode_LBP(int maxIter, double cutoff, int verbose)
{
	double *messages_1 = (double *) R_alloc(maxState * nEdges, sizeof(double)); // Messages from n2 to n1 at edge (n1, n2)
	double *messages_2 = (double *) R_alloc(maxState * nEdges, sizeof(double)); // Messages from n1 to n2 at edge (n1, n2)
	LoopyBP(messages_1, messages_2, maxIter, cutoff, verbose, true);
	Message2NodeBelief(messages_1, messages_2);
	MaxOfMarginals();
}
