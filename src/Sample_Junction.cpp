#include "CRF.h"

SEXP Sample_Junction(SEXP _crf, SEXP _size)
{
	CRF crf(_crf);
	crf.Init_Samples(_size);
	crf.Init_NodeBel();
	crf.Init_EdgeBel();
	crf.Sample_Junction();
	return(crf._samples);
}

void CRF::Sample_Junction(int size)
{
	void *vmax = vmaxget(); 

	JunctionTree tree(*this);
	tree.SendMessages();
	tree.Sample(size);

	vmaxset(vmax); 
}
