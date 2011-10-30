#include "CRF.h"

SEXP Infer_TRBP(SEXP _crf, SEXP _maxIter, SEXP _cutoff, SEXP _verbose)
{
	int maxIter = INTEGER_POINTER(AS_INTEGER(_maxIter))[0];
	double cutoff = NUMERIC_POINTER(AS_NUMERIC(_cutoff))[0];
	int verbose = INTEGER_POINTER(AS_INTEGER(_verbose))[0];

	CRF crf(_crf);
	crf.Init_Belief();
	crf.Infer_TRBP(maxIter, cutoff, verbose);

	return(crf._belief);
}

void CRF::Infer_TRBP(int maxIter, double cutoff, int verbose)
{
	double *messages_1 = (double *) R_alloc(maxState * nEdges, sizeof(double)); // Messages from n2 to n1 at edge (n1, n2)
	double *messages_2 = (double *) R_alloc(maxState * nEdges, sizeof(double)); // Messages from n1 to n2 at edge (n1, n2)
	double *mu = (double *) R_alloc(nEdges, sizeof(double));
	double **scaleEdgePot = (double **) R_alloc(nEdges, sizeof(double *));
	for (int i = 0; i < nEdges; i++)
	{
		scaleEdgePot[i] = (double *) R_alloc(nEdgeStates[i], sizeof(double));
	}
	TRBP_Weights(mu);
	TRBP_ScaleEdgePot(mu, scaleEdgePot);
	TRBP(messages_1, messages_2, mu, scaleEdgePot, maxIter, cutoff, verbose);
	TRBP_Message2NodeBelief(messages_1, messages_2, mu);
	TRBP_Message2EdgeBelief(messages_1, messages_2, mu, scaleEdgePot);
	TRBP_BetheFreeEnergy(mu);
}
