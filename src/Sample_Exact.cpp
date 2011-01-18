#include "CRF.h"

SEXP Sample_Exact(SEXP _crf, SEXP _size)
{
	CRF crf(_crf);
	crf.Init_Samples(_size);
	crf.Sample_Exact();
	return(crf._samples);
}

void CRF::Sample_Exact()
{
	int *y = (int *) R_alloc(nNodes, sizeof(int));
	for (int i = 0; i < nNodes; i++)
		y[i] = 0;

	double pot, Z = 0;
	int index;
	while(1)
	{
		pot = 1;

		/* Node potentials */
		for (int i = 0; i < nNodes; i++)
			pot *= nodePot[i + nNodes * y[i]];

		/* Edge potentials */
		for (int i = 0; i < nEdges; i++)
			pot *= edgePot[y[edges[i]-1] + maxState * (y[edges[i+nEdges]-1] + maxState * i)];

		/* Update Z */
		Z += pot;

		/* Next configuration */
		for (index = 0; index < nNodes; index++)
		{
			y[index] += 1;
			if (y[index] < nStates[index])
				break;
			else
				y[index] = 0;
		}

		if (index == nNodes)
			break;
	}

	/* Sampling */

	double cutoff, cumulativePot;
	GetRNGstate();
	for (int k = 0; k < nSamples; k++)
	{
		for (int i = 0; i < nNodes; i++)
			y[i] = 0;

		cutoff = unif_rand();
		cumulativePot = 0;
		while(1)
		{
			pot = 1;

			/* Node potentials */
			for (int i = 0; i < nNodes; i++)
				pot *= nodePot[i + nNodes * y[i]];

			/* Edge potentials */
			for (int i = 0; i < nEdges; i++)
				pot *= edgePot[y[edges[i]-1] + maxState * (y[edges[i+nEdges]-1] + maxState * i)];

			/* Update cumulative potential */
			cumulativePot += pot;

			if (cumulativePot/Z > cutoff)
				break;

			/* Next configuration */
			for (index = 0; index < nNodes; index++)
			{
				y[index] += 1;
				if (y[index] < nStates[index])
					break;
				else
					y[index] = 0;
			}

			if (index == nNodes)
				break;
		}

		for (int i = 0; i < nNodes; i++)
			samples[k + nSamples * i] = y[i] + 1;
	}
	PutRNGstate();
}
