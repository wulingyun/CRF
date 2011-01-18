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
	int k, maxSample = -1;
	for (int i = 0; i < nSamples; i++)
	{
		pot = 1;
		for (int j = 0; j < nNodes; j++)
		{
			k = j + nNodes * (samples[i + nSamples * j] - 1);
			pot *= nodePot[k];
		}
		for (int j = 0; j < nEdges; j++)
		{
			k = samples[i + nSamples * (edges[j] - 1)] - 1 + maxState * (samples[i + nSamples * (edges[j + nEdges] - 1)] - 1 + maxState * j);
			pot *= edgePot[k];
		}

		if (pot > maxPot)
		{
			maxPot = pot;
			maxSample = i;
		}
	}

	for (int i = 0; i < nNodes; i++)
		labels[i] = samples[maxSample + nSamples * i];
}
