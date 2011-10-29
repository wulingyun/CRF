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
			pot *= EdgePot(j, Samples(i, EdgesBegin(j))-1, Samples(i, EdgesEnd(j))-1);

		if (pot > maxPot)
		{
			maxPot = pot;
			maxSample = i;
		}
	}

	for (int i = 0; i < nNodes; i++)
		labels[i] = Samples(maxSample, i);
}
