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

	MessagesInit();
	TreeBP();
	Messages2NodeBel();
	Messages2EdgeBel();
	BetheFreeEnergy();

	vmaxset(vmax); 
}
