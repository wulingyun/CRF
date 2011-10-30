#include "CRF.h"

SEXP Infer_Cutset(SEXP _crf, SEXP _engine)
{
	CRFclamped crf(_crf);
	crf.Init_Belief();
	crf.original.Init_Belief();
	crf.Infer_Cutset(INTEGER_POINTER(AS_INTEGER(_engine))[0]);
	return(crf.original._belief);
}

void CRFclamped::Infer_Cutset(int engine)
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
	while (1)
	{
		R_CheckUserInterrupt();

		/* Reset node potentials */
		Reset_NodePot();

		/* Infer clamped CRF */
		switch (engine)
		{
		case 0:
			break;
		case 1:
			Infer_Exact();
			break;
		case 2:
			Infer_Chain();
			break;
		case 3:
			Infer_Tree();
			break;
		default:
			Infer_Tree();
			break;
		}

		/* Calculate potential */
		pot = exp(*logZ);
		for (int i = 0; i < original.nNodes; i++)
			if (clamped[i] > 0)
				pot *= original.NodePot(i, y[i]);
		for (int i = 0; i < original.nEdges; i++)
		{
			n1 = original.EdgesBegin(i);
			n2 = original.EdgesEnd(i);
			if (clamped[n1] > 0 && clamped[n2] > 0)
				pot *= original.EdgePot(i, y[n1], y[n2]);
		}

		/* Update node belief */
		for (int i = 0; i < original.nNodes; i++)
		{
			if (clamped[i] > 0)
				original.NodeBel(i, y[i]) += pot;
			else
				for (int j = 0; j < original.nStates[i]; j++)
				{
					original.NodeBel(i, j) += NodeBel(nodeMap[i] - 1, j) * pot;
				}
		}

		/* Update edge belief */
		for (int i = 0; i < original.nEdges; i++)
		{
			n1 = original.EdgesBegin(i);
			n2 = original.EdgesEnd(i);
			if (clamped[n1] > 0)
			{
				if (clamped[n2] > 0)
					original.EdgeBel(i, y[n1], y[n2]) += pot;
				else
					for (int j = 0; j < original.nStates[n2]; j++)
						original.EdgeBel(i, y[n1], j) += NodeBel(nodeMap[n2] - 1, j) * pot;
			}
			else
			{
				if (clamped[n2] > 0)
				{
					for (int j = 0; j < original.nStates[n1]; j++)
						original.EdgeBel(i, j, y[n2]) += NodeBel(nodeMap[n1] - 1, j) * pot;
				}
				else
				{
					for (int j1 = 0; j1 < original.nStates[n1]; j1++)
						for (int j2 = 0; j2 < original.nStates[n2]; j2++)
							original.EdgeBel(i, j1, j2) += EdgeBel(edgeMap[i] - 1, j1, j2) * pot;
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
	for (int i = 0; i < original.nEdges; i++)
	{
		for (int j = 0; j < original.nEdgeStates[i]; j++)
			original.edgeBel[i][j] /= Z;
	}
	*(original.logZ) = log(Z);
}
