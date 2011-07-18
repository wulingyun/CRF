#include "CRF.h"

SEXP Decode_ICM(SEXP _crf, SEXP _start, SEXP _restart, SEXP _maxIter)
{
	PROTECT(_start = AS_INTEGER(_start));
	int *start = INTEGER_POINTER(_start);
	bool restart = LOGICAL_POINTER(AS_LOGICAL(_restart))[0];
	int maxIter = INTEGER_POINTER(AS_INTEGER(_maxIter))[0];

	CRF crf(_crf);
	crf.Init_Labels();
	crf.Decode_ICM(start, restart, maxIter);

	UNPROTECT_PTR(_start);
	return(crf._labels);
}

void CRF::Decode_ICM(int *start, bool restart, int maxIter)
{
	int *y = (int *) R_alloc(nNodes, sizeof(int));
	double max;
	if (start)
		for (int i = 0; i < nNodes; i++)
			y[i] = start[i] - 1;
	else
		for (int i = 0; i < nNodes; i++)
		{
			max = -1;
			for (int j = 0; j < nStates[i]; j++)
				if (max < NodePot(i,j))
				{
					max = NodePot(i,j);
					y[i] = j;
				}
		}

	double *pot = (double *) R_alloc(maxState, sizeof(double));

	int e, n1, n2;
	bool done;

	double Z, maxZ = Get_Potential(y);
	for (int i = 0; i < nNodes; i++)
		labels[i] = y[i] + 1;

	GetRNGstate();
	for (int iter = 0; iter < maxIter; iter++)
	{
		done = false;
		while (!done)
		{
			R_CheckUserInterrupt();

			done = true;
			for (int i = 0; i < nNodes; i++)
			{
				for (int j = 0; j < nStates[i]; j++)
					pot[j] = NodePot(i,j);

				for (int j = 0; j < nAdj[i]; j++)
				{
					e = adjEdges[i][j] - 1;
					n1 = Edges(e,0) - 1;
					n2 = Edges(e,1) - 1;
					if (i == n1)
						for (int k = 0; k < nStates[i]; k++)
							pot[k] *= EdgePot(k, y[n2], e);
					else
						for (int k = 0; k < nStates[i]; k++)
							pot[k] *= EdgePot(y[n1], k, e);
				}

				max = pot[y[i]];
				for (int j = 0; j < nStates[i]; j++)
				{
					if (max < pot[j])
					{
						max = pot[j];
						y[i] = j;
						done = 0;
					}
				}
			}
		}

		Z = Get_Potential(y);
		if (Z > maxZ)
		{
			maxZ = Z;
			for (int i = 0; i < nNodes; i++)
				labels[i] = y[i] + 1;
		}

		if (restart)
			for (int i = 0; i < nNodes; i++)
				y[i] = ceil(unif_rand() * nStates[i]) - 1;
		else
			break;
	}
	PutRNGstate();
}
