#include "CRF.h"

SEXP Sample_Gibbs(SEXP _crf, SEXP _size, SEXP _burnIn, SEXP _start)
{
	int burnIn = INTEGER_POINTER(AS_INTEGER(_burnIn))[0];
	PROTECT(_start = AS_INTEGER(_start));
	int *start = INTEGER_POINTER(_start);

	CRF crf(_crf);
	crf.Init_Samples(_size);
	crf.Sample_Gibbs(burnIn, start);

	UNPROTECT_PTR(_start);
	return(crf._samples);
}

void CRF::Sample_Gibbs(int burnIn, int *start, int size)
{
	if (size <= 0)
		size = nSamples;
	else if (size > nSamples)
		Init_Samples(size);

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

	double sumProb, *prob = (double *) R_alloc(maxState, sizeof(double));

	int e, n, n1, n2;
	double *p_nodePot, *p_edgePot;

	GetRNGstate();
	for (int iter = 0; iter < burnIn+size; iter++)
	{
		R_CheckUserInterrupt();

		for (int i = 0; i < nNodes; i++)
		{
			n = nStates[i];
			p_nodePot = nodePot + i;
			for (int j = 0; j < n; j++)
			{
				prob[j] = p_nodePot[0];
				p_nodePot += nNodes;
			}
			for (int j = 0; j < nAdj[i]; j++)
			{
				e = adjEdges[i][j] - 1;
				n1 = edges[e] - 1;
				n2 = edges[e + nEdges] - 1;
				if (n1 == i)
				{
					p_edgePot = edgePot + maxState * (y[n2] + maxState * e);
					for (int k = 0; k < n; k++)
						prob[k] *= p_edgePot[k];
				}
				else
				{
					p_edgePot = edgePot + y[n1] + maxState * maxState * e;
					for (int k = 0; k < n; k++)
					{
						prob[k] *= p_edgePot[0];
						p_edgePot += maxState;
					}
				}
			}
			sumProb = 0;
			for (int j = 0; j < n; j++)
				sumProb += prob[j];
			for (int j = 0; j < n; j++)
				prob[j] /= sumProb;

			y[i] = sample(n, prob);
		}

		if (iter >= burnIn)
			for (int i = 0; i < nNodes; i++)
				Samples(iter - burnIn, i) = y[i] + 1;
	}
	PutRNGstate();
}
