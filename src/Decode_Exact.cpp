#include "CRF.h"

SEXP Decode_Exact(SEXP _crf)
{
	CRF crf(_crf);
	crf.Init_Labels();
	crf.Decode_Exact();
	return(crf._labels);
}

void CRF::Decode_Exact()
{
	int *y = (int *) R_alloc(nNodes, sizeof(int));
	for (int i = 0; i < nNodes; i++)
		y[i] = 0;

	double pot, maxPot = -1;
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

		/* Record the best potentials */
		if (pot > maxPot)
		{
			maxPot = pot;
			for (int i = 0; i < nNodes; i++)
				labels[i] = y[i] + 1;
		}

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
}
