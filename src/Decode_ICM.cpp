#include "CRF.h"

SEXP Decode_ICM(SEXP _crf, SEXP _restart, SEXP _start)
{
	int restart = INTEGER_POINTER(AS_INTEGER(_restart))[0];
	PROTECT(_start = AS_INTEGER(_start));
	int *start = INTEGER_POINTER(_start);

	CRF crf(_crf);
	crf.Init_Labels();
	crf.Decode_ICM(restart, start);

	UNPROTECT_PTR(_start);
	return(crf._labels);
}

void CRF::Decode_ICM(int restart, int *start)
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

	int e, n1, n2;
	bool done;

	GetRNGstate();
	for (int iter = 0; iter <= restart; iter++)
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
					e = AdjEdges(i, j);
					n1 = EdgesBegin(e);
					n2 = EdgesEnd(e);
					if (i == n1)
						for (int k = 0; k < nStates[i]; k++)
							pot[k] *= EdgePot(e, k, y[n2]);
					else
						for (int k = 0; k < nStates[i]; k++)
							pot[k] *= EdgePot(e, y[n1], k);
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

		if (iter < restart)
			for (int i = 0; i < nNodes; i++)
				y[i] = ceil(unif_rand() * nStates[i]) - 1;
	}
	PutRNGstate();
}
