#include "CRF.h"

SEXP Infer_LBP(SEXP _crf, SEXP _maxIter, SEXP _cutoff, SEXP _verbose)
{
	PROTECT(_maxIter = AS_INTEGER(_maxIter));
	int maxIter = INTEGER_POINTER(_maxIter)[0];
	PROTECT(_cutoff = AS_NUMERIC(_cutoff));
	double cutoff = NUMERIC_POINTER(_cutoff)[0];
	PROTECT(_verbose = AS_INTEGER(_verbose));
	int verbose = INTEGER_POINTER(_verbose)[0];

	CRF crf(_crf);
	crf.Init_Belief();
	crf.Infer_LBP(maxIter, cutoff, verbose);

	UNPROTECT(3);
	return(crf._belief);
}

void CRF::Infer_LBP(int maxIter, double cutoff, int verbose)
{
	double *messages_1 = (double *) R_alloc(maxState * nEdges, sizeof(double)); // Messages from n2 to n1 at edge (n1, n2)
	double *messages_2 = (double *) R_alloc(maxState * nEdges, sizeof(double)); // Messages from n1 to n2 at edge (n1, n2)
	LoopyBP(messages_1, messages_2, maxIter, cutoff, verbose);
	Message2NodeBelief(messages_1, messages_2);
	Message2EdgeBelief(messages_1, messages_2);
	BetheFreeEnergy();
}
