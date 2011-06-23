#include "CRF.h"

SEXP Infer_Tree(SEXP _crf)
{
	CRF crf(_crf);
	crf.Init_Belief();
	crf.Infer_Tree();
	return(crf._belief);
}

void CRF::Infer_Tree()
{
	void *vmax = vmaxget(); 

	double *messages_1 = (double *) R_alloc(maxState * nEdges, sizeof(double)); // Messages from n2 to n1 at edge (n1, n2)
	double *messages_2 = (double *) R_alloc(maxState * nEdges, sizeof(double)); // Messages from n1 to n2 at edge (n1, n2)
	TreeBP(messages_1, messages_2);
	Message2NodeBelief(messages_1, messages_2);
	Message2EdgeBelief(messages_1, messages_2);
	BetheFreeEnergy();

	vmaxset(vmax); 
}
