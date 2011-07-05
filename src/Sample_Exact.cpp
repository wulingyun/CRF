#include "CRF.h"

SEXP Sample_Exact(SEXP _crf, SEXP _size)
{
	CRF crf(_crf);
	crf.Init_Samples(_size);
	crf.Sample_Exact();
	return(crf._samples);
}

void CRF::Sample_Exact(int size)
{
	if (size <= 0)
		size = nSamples;
	else if (size > nSamples)
		Init_Samples(size);

	int *y = (int *) R_alloc(nNodes, sizeof(int));
	for (int i = 0; i < nNodes; i++)
		y[i] = 0;

	double Z = 0;
	int index;
	while (1)
	{
		R_CheckUserInterrupt();

		/* Calculate potential and update Z */
		Z += Get_Potential(y);

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

	double *cutoff = (double *) R_alloc(size, sizeof(double));
	GetRNGstate();
	for (int k = 0; k < size; k++)
		cutoff[k] = unif_rand() * Z;
	PutRNGstate();

	for (int i = 0; i < nNodes; i++)
		y[i] = 0;

	int remain = size;
	double done = Z * 10;
	double cumulativePot = 0;
	while (1)
	{
		R_CheckUserInterrupt();

		/* Update cumulative potential */
		cumulativePot += Get_Potential(y);

		for (int k = 0; k < size; k++)
			if (cumulativePot > cutoff[k])
			{
				for (int i = 0; i < nNodes; i++)
					Samples(k, i) = y[i] + 1;
				cutoff[k] = done;
				remain--;
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

		if (index == nNodes || remain <= 0)
			break;
	}
}
