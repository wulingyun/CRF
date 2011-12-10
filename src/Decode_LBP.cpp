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
	MessagesInit();
	LoopyBP(maxIter, cutoff, verbose, true);
	Messages2NodeBel();
	MaxOfMarginals();
}
