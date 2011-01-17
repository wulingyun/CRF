#include "CRF.h"

SEXP Decode_Tree(SEXP _crf)
{
	CRF crf(_crf);

	SEXP _labels;
	PROTECT(_labels = NEW_INTEGER(crf.nNodes));
	int *labels = INTEGER_POINTER(_labels);
	setValues(_labels, labels, -1);

	crf.Decode_Tree(labels);

	UNPROTECT(1);

	return(_labels);
}

void CRF::Decode_Tree(int *labels)
{
	double *messages_1 = (double *) R_alloc(maxState * nEdges, sizeof(double)); // Messages from n2 to n1 at edge (n1, n2)
	double *messages_2 = (double *) R_alloc(maxState * nEdges, sizeof(double)); // Messages from n1 to n2 at edge (n1, n2)
	double *nodeBel = (double *) R_alloc(nNodes * maxState, sizeof(double));
	TreeBP_max(messages_1, messages_2);
	Message2NodeBelief(messages_1, messages_2, nodeBel);
	MaxOfMarginals(nodeBel, labels);
}
