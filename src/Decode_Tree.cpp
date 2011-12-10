#include "CRF.h"

SEXP Decode_Tree(SEXP _crf)
{
	CRF crf(_crf);
	crf.Init_Labels();
	crf.Init_NodeBel();
	crf.Decode_Tree();
	return(crf._labels);
}

void CRF::Decode_Tree()
{
	void *vmax = vmaxget(); 

	MessagesInit();
	TreeBP(true);
	Messages2NodeBel();
	MaxOfMarginals();

	vmaxset(vmax); 
}
