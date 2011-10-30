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
	int k, s1, s2, maxSample = -1;
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
			s1 = Samples(i, EdgesBegin(j)) - 1;
			s2 = Samples(i, EdgesEnd(j)) - 1;
			EdgeBel(j, s1, s2)++;
			pot *= EdgePot(j, s1, s2);
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
	for (int i = 0; i < nEdges; i++)
	{
		for (int j = 0; j < nEdgeStates[i]; j++)
			edgeBel[i][j] /= nSamples;
	}
	*logZ = log(maxPot * nSamples / freq);
}
