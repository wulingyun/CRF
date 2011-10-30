#include "CRF.h"

SEXP Infer_Exact(SEXP _crf)
{
	CRF crf(_crf);
	crf.Init_Belief();
	crf.Infer_Exact();
	return(crf._belief);
}

void CRF::Infer_Exact()
{
	int *y = (int *) R_alloc(nNodes, sizeof(int));
	for (int i = 0; i < nNodes; i++)
		y[i] = 0;

	double pot, Z = 0;
	int index;
	while (1)
	{
		R_CheckUserInterrupt();

		/* Calculate potential */
		pot = Get_Potential(y);

		/* Node belief */
		for (int i = 0; i < nNodes; i++)
			nodeBel[i + nNodes * y[i]] += pot;

		/* Edge belief */
		for (int i = 0; i < nEdges; i++)
			EdgeBel(i, y[EdgesBegin(i)], y[EdgesEnd(i)]) += pot;

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

	/* Normalization */
	for (int i = 0; i < length(_nodeBel); i++)
		nodeBel[i] /= Z;
	for (int i = 0; i < nEdges; i++)
	{
		for (int j = 0; j < nEdgeStates[i]; j++)
			edgeBel[i][j] /= Z;
	}
	*logZ = log(Z);
}
