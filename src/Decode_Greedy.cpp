#include "CRF.h"

SEXP Decode_Greedy(SEXP _crf, SEXP _restart, SEXP _start)
{
	int restart = INTEGER_POINTER(AS_INTEGER(_restart))[0];
	PROTECT(_start = AS_INTEGER(_start));
	int *start = INTEGER_POINTER(_start);

	CRF crf(_crf);
	crf.Init_Labels();
	crf.Decode_Greedy(restart, start);

	UNPROTECT_PTR(_start);
	return(crf._labels);
}

void CRF::Decode_Greedy(int restart, int *start)
{
	if (restart < 0)
		restart = 0;

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

	double Z, maxZ = Get_Potential(y);
	for (int i = 0; i < nNodes; i++)
		labels[i] = y[i] + 1;

	double *pot = (double *) R_alloc(maxState, sizeof(double));
	double *diff = (double *) R_alloc(nNodes, sizeof(double));
	int *move = (int *) R_alloc(nNodes, sizeof(double));

	double ref;
	int e, n1, n2, index;

	GetRNGstate();
	for (int iter = 0; iter <= restart; iter++)
	{
		while (1)
		{
			R_CheckUserInterrupt();

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

				ref = pot[y[i]];
				if (ref != 0)
					for (int j = 0; j < nStates[i]; j++)
						pot[j] /= ref;

				diff[i] = -1;
				for (int j = 0; j < nStates[i]; j++)
				{
					if (diff[i] < pot[j])
					{
						diff[i] = pot[j];
						move[i] = j;
					}
				}
			}

			max = -1;
			index = -1;
			for (int i = 0; i < nNodes; i++)
			{
				if (max < diff[i])
				{
					max = diff[i];
					index = i;
				}
			}

			if (max <= 1)
				break;
			else
				y[index] = move[index];
		}

		Z = Get_Potential(y);
		if (Z > maxZ)
		{
			maxZ = Z;
			for (int i = 0; i < nNodes; i++)
				labels[i] = y[i] + 1;
		}

		if (iter < restart)
			for (int i = 0; i < nNodes; i++)
				y[i] = ceil(unif_rand() * nStates[i]) - 1;
	}
	PutRNGstate();
}
