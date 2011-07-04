#include "CRF.h"

SEXP Infer_Sample(SEXP _crf, SEXP _samples)
{
	CRF crf(_crf);
	crf.Init_Belief();
	crf.Set_Samples(_samples);
	crf.Infer_Sample();
	return(crf._belief);
}

void CRF::Infer_Sample()
{
	double pot, maxPot = -1;
	int k, maxSample = -1;
	for (int i = 0; i < nSamples; i++)
	{
		R_CheckUserInterrupt();

		pot = 1;
		for (int j = 0; j < nNodes; j++)
		{
			k = j + nNodes * (Samples(i, j) - 1);
			nodeBel[k]++;
			pot *= nodePot[k];
		}
		for (int j = 0; j < nEdges; j++)
		{
			k = Samples(i, Edges(j,0)-1) - 1 + maxState * (Samples(i, Edges(j,1)-1) - 1 + maxState * j);
			edgeBel[k]++;
			pot *= edgePot[k];
		}

		if (pot > maxPot)
		{
			maxPot = pot;
			maxSample = i;
		}
	}
	int same, freq = 0;
	for (int i = 0; i < nSamples; i++)
	{
		R_CheckUserInterrupt();

		same = 1;
		for (int j = 0; j < nNodes; j++)
		{
			if (Samples(i, j) != Samples(maxSample, j))
			{
				same = 0;
				break;
			}
		}
		if (same)
			freq++;
	}
	for (int i = 0; i < length(_nodeBel); i++)
		nodeBel[i] /= nSamples;
	for (int i = 0; i < length(_edgeBel); i++)
		edgeBel[i] /= nSamples;
	*logZ = log(maxPot * nSamples / freq);
}
