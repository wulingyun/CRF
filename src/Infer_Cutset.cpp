#include "CRF.h"

SEXP Infer_Cutset(SEXP _crf)
{
	CRFclamped crf(_crf);
	crf.Init_Belief();
	crf.original.Init_Belief();
	crf.Infer_Cutset();
	return(crf.original._belief);
}

void CRFclamped::Infer_Cutset()
{
	int *y = (int *) R_alloc(original.nNodes, sizeof(int));
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

	double pot, Z = 0;
	int index, n1, n2;
	while(1)
	{
		/* Reset node potentials */
		Reset_NodePot();

		/* Infer clamped CRF */
		Infer_Tree();

		/* Calculate potential */
		pot = exp(*logZ);
		for (int i = 0; i < original.nNodes; i++)
			if (clamped[i] > 0)
				pot *= original.nodePot[i + original.nNodes * y[i]];
		for (int i = 0; i < original.nEdges; i++)
		{
			n1 = original.edges[i] - 1;
			n2 = original.edges[i + original.nEdges] - 1;
			if (clamped[n1] > 0 && clamped[n2] > 0)
				pot *= original.edgePot[y[n1] + original.maxState * (y[n2] + original.maxState * i)];
		}

		/* Update node belief */
		for (int i = 0; i < original.nNodes; i++)
		{
			if (clamped[i] > 0)
				original.nodeBel[i + original.nNodes * y[i]] += pot;
			else
				for (int j = 0; j < original.maxState; j++)
				{
					original.nodeBel[i + original.nNodes * j] += nodeBel[nodeMap[i] - 1 + nNodes * j] * pot;
				}
		}

		/* Update edge belief */
		for (int i = 0; i < original.nEdges; i++)
		{
			n1 = original.edges[i] - 1;
			n2 = original.edges[i + original.nEdges] - 1;
			if (clamped[n1] > 0)
			{
				if (clamped[n2] > 0)
					original.edgeBel[y[n1] + original.maxState * (y[n2] + original.maxState * i)] += pot;
				else
					for (int j = 0; j < original.maxState; j++)
					{
						original.edgeBel[y[n1] + original.maxState * (j + original.maxState * i)] += nodeBel[nodeMap[n2] - 1 + nNodes * j] * pot;
					}
			}
			else
			{
				if (clamped[n2] > 0)
				{
					for (int j = 0; j < original.maxState; j++)
					{
						original.edgeBel[j + original.maxState * (y[n2] + original.maxState * i)] += nodeBel[nodeMap[n1] - 1 + nNodes * j] * pot;
					}
				}
				else
				{
					for (int j1 = 0; j1 < original.maxState; j1++)
						for (int j2 = 0; j2 < original.maxState; j2++)
						{
							original.edgeBel[j1 + original.maxState * (j2 + original.maxState * i)] += edgeBel[j1 + maxState * (j2 + maxState * (edgeMap[i] - 1))] * pot;
						}
				}
			}
		}

		/* Update Z */
		Z += pot;

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

	/* Normalization */
	for (int i = 0; i < length(original._nodeBel); i++)
		original.nodeBel[i] /= Z;
	for (int i = 0; i < length(original._edgeBel); i++)
		original.edgeBel[i] /= Z;
	*(original.logZ) = log(Z);
}
