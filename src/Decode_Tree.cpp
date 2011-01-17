#include "CRF.h"

SEXP Decode_Tree(SEXP _crf)
{
	CRF crf(_crf);
	crf.Init_Decoding2();
	crf.Decode_Tree();
	return(crf._labels);
}

void CRF::Decode_Tree()
{
	double *messages_1 = (double *) R_alloc(maxState * nEdges, sizeof(double)); // Messages from n2 to n1 at edge (n1, n2)
	double *messages_2 = (double *) R_alloc(maxState * nEdges, sizeof(double)); // Messages from n1 to n2 at edge (n1, n2)
	TreeBP_max(messages_1, messages_2);
	Message2NodeBelief(messages_1, messages_2);
	MaxOfMarginals();
}
