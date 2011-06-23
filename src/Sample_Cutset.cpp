#include "CRF.h"

SEXP Sample_Cutset(SEXP _crf, SEXP _size)
{
	CRFclamped crf(_crf);
	crf.Init_Belief();
	crf.Init_Samples(1);
	crf.original.Init_Samples(_size);
	crf.Sample_Cutset();
	return(crf.original._samples);
}

void CRFclamped::Sample_Cutset()
{
	int *y = (int *) R_alloc(original.nNodes, sizeof(int));
	int nPot = 1;
	for (int i = 0; i < original.nNodes; i++)
	{
		if (clamped[i] > 0)
		{
			nPot *= original.nStates[i];
			clamped[i] = 1;
			y[i] = 0;
		}
		else
		{
			clamped[i] = 0;
			y[i] = -1;
		}
	}

	double *pot = (double *) R_alloc(nPot, sizeof(double));
	double Z = 0;
	int n = 0;
	int index, n1, n2;
	while(1)
	{
		/* Reset node potentials */
		Reset_NodePot();

		/* Infer clamped CRF */
		Infer_Tree();

		/* Calculate potential */
		pot[n] = exp(*logZ);
		for (int i = 0; i < original.nNodes; i++)
			if (clamped[i] > 0)
				pot[n] *= original.NodePot(i, y[i]);
		for (int i = 0; i < original.nEdges; i++)
		{
			n1 = original.edges[i] - 1;
			n2 = original.edges[i + original.nEdges] - 1;
			if (clamped[n1] > 0 && clamped[n2] > 0)
				pot[n] *= original.EdgePot(y[n1], y[n2], i);
		}

		/* Update Z */
		Z += pot[n++];

		/* Next configuration */
		for (index = 0; index < original.nNodes; index++)
		{
			if (clamped[index] == 0)
				continue;
			clamped[index]++;
			y[index]++;
			if (y[index] < original.nStates[index])
				break;
			else
			{
				clamped[index] = 1;
				y[index] = 0;
			}
		}

		if (index == original.nNodes)
			break;
	}

	/* Sampling */

	double *cutoff = (double *) R_alloc(original.nSamples, sizeof(double));
	GetRNGstate();
	for (int k = 0; k < original.nSamples; k++)
		cutoff[k] = unif_rand() * Z;
	PutRNGstate();

	for (int i = 0; i < original.nNodes; i++)
	{
		if (clamped[i] > 0)
		{
			clamped[i] = 1;
			y[i] = 0;
		}
		else
		{
			clamped[i] = 0;
			y[i] = -1;
		}
	}

	int remain = original.nSamples;
	double done = Z * 10;
	double cumulativePot = 0;
	n = 0;
	int m;
	while(1)
	{
		/* Reset node potentials */
		Reset_NodePot();

		/* Update cumulative potential */
		cumulativePot += pot[n++];

		/* Count samples */
		m = 0;
		for (int k = 0; k < original.nSamples; k++)
			if (cumulativePot > cutoff[k])
				m++;

		/* Generate samples */
		if (m > 0)
		{
			Sample_Tree(m);
			m = 0;
			for (int k = 0; k < original.nSamples; k++)
				if (cumulativePot > cutoff[k])
				{
					for (int i = 0; i < original.nNodes; i++)
					{
						if (clamped[i] > 0)
							original.Samples(k, i) = clamped[i];
						else
							original.Samples(k, i) = Samples(m, nodeMap[i] - 1);
					}
					cutoff[k] = done;
					remain--;
					m++;
				}
		}

		/* Next configuration */
		for (index = 0; index < original.nNodes; index++)
		{
			if (clamped[index] == 0)
				continue;
			clamped[index]++;
			y[index]++;
			if (y[index] < original.nStates[index])
				break;
			else
			{
				clamped[index] = 1;
				y[index] = 0;
			}
		}

		if (index == original.nNodes || remain <= 0)
			break;
	}
}
