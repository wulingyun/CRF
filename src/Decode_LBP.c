#include "CRF.h"

SEXP Decode_LBP(SEXP _crf, SEXP _maxIter, SEXP _cutoff, SEXP _verbose)
{
	CRFinfo crf;
	openCRF(&crf, _crf);

	PROTECT(_maxIter = AS_INTEGER(_maxIter));
	int maxIter = INTEGER_POINTER(_maxIter)[0];
	PROTECT(_cutoff = AS_NUMERIC(_cutoff));
	double cutoff = NUMERIC_POINTER(_cutoff)[0];
	PROTECT(_verbose = AS_INTEGER(_verbose));
	int verbose = INTEGER_POINTER(_verbose)[0];

	SEXP _labels;
	PROTECT(_labels = NEW_INTEGER(crf.nNodes));
	int *labels = INTEGER_POINTER(_labels);
	setValues(_labels, labels, -1);

	_Decode_LBP(&crf, labels, maxIter, cutoff, verbose);

	UNPROTECT(4);
	closeCRF(&crf);

	return(_labels);
}

void _Decode_LBP(CRFinfo *crf, int *labels, int maxIter, double cutoff, int verbose)
{
	/* Loopy BP */

	double *messages_1 = (double *) R_alloc(crf->maxState * crf->nEdges, sizeof(double)); // Messages from n2 to n1 at edge (n1, n2)
	double *messages_2 = (double *) R_alloc(crf->maxState * crf->nEdges, sizeof(double)); // Messages from n1 to n2 at edge (n1, n2)
	LoopyBP_max(crf, messages_1, messages_2, maxIter, cutoff, verbose);

	/* Node beliefs */

	double *nodeBel = (double *) R_alloc(crf->nNodes * crf->maxState, sizeof(double));
	Message2NodeBelief(crf, messages_1, messages_2, nodeBel);

	/* Labels */

	double maxBel, *p_nodeBel;
	for (int i = 0; i < crf->nNodes; i++)
	{
		maxBel = -1;
		p_nodeBel = nodeBel + i;
		for (int j = 0; j < crf->nStates[i]; j++)
		{
			if (p_nodeBel[0] > maxBel)
			{
				maxBel = p_nodeBel[0];
				labels[i] = j;
			}
			p_nodeBel += crf->nNodes;
		}
	}

	for (int i = 0; i < crf->nNodes; i++)
		labels[i]++;
}
