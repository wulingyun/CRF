#include "CRF.h"

SEXP Decode_Junction(SEXP _crf)
{
	CRF crf(_crf);
	crf.Init_Labels();
	crf.Init_NodeBel();
	crf.Decode_Junction();
	return(crf._labels);
}

void CRF::Decode_Junction()
{
	void *vmax = vmaxget(); 

	jTree = new JunctionTree(*this);

	vmaxset(vmax); 
}
