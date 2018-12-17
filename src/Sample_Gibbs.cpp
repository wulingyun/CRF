#include "CRF.h"

SEXP Sample_Gibbs(SEXP _crf, SEXP _size, SEXP _burnIn, SEXP _start)
{
	int burnIn = INTEGER_POINTER(AS_INTEGER(_burnIn))[0];

	CRF crf(_crf);
	crf.Init_Samples(_size);

	PROTECT(_start = AS_INTEGER(_start));
	int *start = INTEGER_POINTER(_start);

	crf.Sample_Gibbs(burnIn, start);

	UNPROTECT(1);

	return(crf._samples);
}

void CRF::Sample_Gibbs(int burnIn, int *start)
{
	int	size = nSamples;

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

	GetRNGstate();
	for (int iter = 0; iter < burnIn+size; iter++)
	{
		R_CheckUserInterrupt();

		for (int i = 0; i < nNodes; i++)
		{
			n = nStates[i];
			for (int j = 0; j < n; j++)
				prob[j] = NodePot(i, j);
			for (int j = 0; j < nAdj[i]; j++)
			{
				e = AdjEdges(i, j);
				n1 = EdgesBegin(e);
				n2 = EdgesEnd(e);
				if (n1 == i)
				{
					for (int k = 0; k < n; k++)
						prob[k] *= EdgePot(e, k, y[n2]);
				}
				else
				{
					for (int k = 0; k < n; k++)
						prob[k] *= EdgePot(e, y[n1], k);
				}
			}
			sumProb = 0;
			for (int j = 0; j < n; j++)
				sumProb += prob[j];
			for (int j = 0; j < n; j++)
				prob[j] /= sumProb;

			y[i] = SampleFrom(n, prob);
		}

		if (iter >= burnIn)
			for (int i = 0; i < nNodes; i++)
				Samples(iter - burnIn, i) = y[i] + 1;
	}
	PutRNGstate();
}
