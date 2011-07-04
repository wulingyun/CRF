#include "CRF.h"

SEXP Decode_Sample(SEXP _crf, SEXP _samples)
{
	CRF crf(_crf);
	crf.Init_Labels();
	crf.Set_Samples(_samples);
	crf.Decode_Sample();
	return(crf._labels);
}

void CRF::Decode_Sample()
{
	double pot, maxPot = -1;
	int maxSample = -1;
	for (int i = 0; i < nSamples; i++)
	{
		R_CheckUserInterrupt();

		pot = 1;
		for (int j = 0; j < nNodes; j++)
			pot *= NodePot(j, Samples(i, j) - 1);
		for (int j = 0; j < nEdges; j++)
			pot *= EdgePot(Samples(i, Edges(j,0)-1) - 1, Samples(i, Edges(j,1)-1) - 1, j);

		if (pot > maxPot)
		{
			maxPot = pot;
			maxSample = i;
		}
	}

	for (int i = 0; i < nNodes; i++)
		labels[i] = Samples(maxSample, i);
}
