#include "CRF.h"

SEXP Infer_Junction(SEXP _crf)
{
	CRF crf(_crf);
	crf.Init_Belief();
	crf.Infer_Junction();
	return(crf._belief);
}

void CRF::Infer_Junction()
{
	void *vmax = vmaxget(); 

	JunctionTree tree(*this);
	tree.SendMessages();
	tree.Messages2EdgeBel();
	BetheFreeEnergy();

	vmaxset(vmax); 
}
