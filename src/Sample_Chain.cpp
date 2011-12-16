#include "CRF.h"

SEXP Sample_Chain(SEXP _crf, SEXP _size)
{
	CRF crf(_crf);
	crf.Init_Samples(_size);
	crf.Sample_Chain();
	return(crf._samples);
}

void CRF::Sample_Chain(int size)
{
	if (size <= 0)
		size = nSamples;
	else if (size > nSamples)
		Init_Samples(size);

	int *y = (int *) R_alloc(nNodes, sizeof(int));
	for (int i = 0; i < nNodes; i++)
		y[i] = 0;

	/* forward pass */

	double *alpha = (double *) R_alloc(nNodes * maxState, sizeof(double));
	for (int i = 0; i < nNodes * maxState; i++)
		alpha[i] = 0;
	double *kappa = (double *) R_alloc(nNodes, sizeof(double));
	for (int i = 0; i < nNodes; i++)
		kappa[i] = 0;

	for (int i = 0; i < nStates[0]; i++)
	{
		alpha[nNodes * i] = nodePot[nNodes * i];
		kappa[0] += alpha[nNodes * i];
	}
	if (kappa[0] != 0)
	{
		for (int i = 0; i < nStates[0]; i++)
			alpha[nNodes * i] /= kappa[0];
	}

	double *p_alpha, *p0_alpha;
	double sumPot;
	for (int i = 1; i < nNodes; i++)
	{
		p0_alpha = alpha + i;
		for (int j = 0; j < nStates[i]; j++)
		{
			p_alpha = alpha + i - 1;
			sumPot = 0;
			for (int k=0; k < nStates[i-1]; k++)
			{
				sumPot += p_alpha[0] * EdgePot(i-1, k, j);
				p_alpha += nNodes;
			}
			p0_alpha[0] = sumPot * NodePot(i, j);
			kappa[i] += p0_alpha[0];
			p0_alpha += nNodes;
		}
		if (kappa[i] != 0)
		{
			p_alpha = alpha + i;
			for (int j = 0; j < nStates[i]; j++)
			{
				p_alpha[0] /= kappa[i];
				p_alpha += nNodes;
			}
		}
	}

	/* backward pass */

	double sumProb, *prob = (double *) R_alloc(maxState, sizeof(double));

	GetRNGstate();
	for (int i = 0; i < size; i++)
	{
		p_alpha = alpha + nNodes - 1;
		for (int j = 0; j < nStates[nNodes-1]; j++)
		{
			prob[j] = p_alpha[0];
			p_alpha += nNodes;
		}
		y[nNodes-1] = SampleFrom(nStates[nNodes-1], prob);
		for (int j = nNodes-2; j >= 0; j--)
		{
			p_alpha = alpha + j;
			sumProb = 0;
			for (int k = 0; k < nStates[j]; k++)
			{
				prob[k] = p_alpha[0] * EdgePot(j, k, y[j+1]);
				sumProb += prob[k];
				p_alpha += nNodes;
			}
			if (sumProb != 0)
			{
				for (int k = 0; k < nStates[j]; k++)
					prob[k] /= sumProb;
			}
			y[j] = SampleFrom(nStates[j], prob);
		}

		for (int j = 0; j < nNodes; j++)
			Samples(i, j) = y[j] + 1;
	}
	PutRNGstate();
}
