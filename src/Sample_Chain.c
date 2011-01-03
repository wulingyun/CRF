#include <time.h>
#include "CRF.h"

SEXP Sample_Chain(SEXP _crf, SEXP _size)
{
	SEXP _nNodes, _nStates, _maxState;
	PROTECT(_nNodes = AS_INTEGER(getListElement(_crf, "n.nodes")));
	PROTECT(_nStates = AS_INTEGER(getListElement(_crf, "n.states")));
	PROTECT(_maxState = AS_INTEGER(getListElement(_crf, "max.state")));
	int nNodes = INTEGER_POINTER(_nNodes)[0];
	int *nStates = INTEGER_POINTER(_nStates);
	int maxState = INTEGER_POINTER(_maxState)[0];

	SEXP _nodePot, _edgePot;
	PROTECT(_nodePot = AS_NUMERIC(getListElement(_crf, "node.pot")));
	PROTECT(_edgePot = AS_NUMERIC(getListElement(_crf, "edge.pot")));
	double *nodePot = NUMERIC_POINTER(_nodePot);
	double *edgePot = NUMERIC_POINTER(_edgePot);

	PROTECT(_size = AS_INTEGER(_size));
	int size = INTEGER_POINTER(_size)[0];

	SEXP _samples;
	PROTECT(_samples = NEW_INTEGER(size * nNodes));
	setDim2(_samples, size, nNodes);
	int *samples = INTEGER_POINTER(_samples);
	setValues(_samples, samples, 0);

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
	for (int i = 0; i < nStates[0]; i++)
		alpha[nNodes * i] /= kappa[0];

	double *p_alpha, *p0_alpha, *p_nodePot, *p_edgePot;
	double sumPot;
	for (int i = 1; i < nNodes; i++)
	{
		p0_alpha = alpha + i;
		p_nodePot = nodePot + i;
		p_edgePot = edgePot + maxState * maxState * (i-1);
		for (int j = 0; j < nStates[i]; j++)
		{
			p_alpha = alpha + i - 1;
			sumPot = 0;
			for (int k=0; k < nStates[i-1]; k++)
			{
				sumPot += p_alpha[0] * p_edgePot[k];
				p_alpha += nNodes;
			}
			p0_alpha[0] = sumPot * p_nodePot[0];
			kappa[i] += p0_alpha[0];
			p0_alpha += nNodes;
			p_nodePot += nNodes;
			p_edgePot += maxState;
		}
		p_alpha = alpha + i;
		for (int j = 0; j < nStates[i]; j++)
		{
			p_alpha[0] /= kappa[i];
			p_alpha += nNodes;
		}
	}

	/* backward pass */

	double *prob = (double *) R_alloc(maxState, sizeof(double));

	srand((int) time(0));
	for (int i = 0; i < size; i++)
	{
		p_alpha = alpha + nNodes - 1;
		for (int j = 0; j < nStates[nNodes-1]; j++)
		{
			prob[j] = p_alpha[0];
			p_alpha += nNodes;
		}
		y[nNodes-1] = sample(nStates[nNodes-1], prob);
		for (int j = nNodes-2; j >= 0; j--)
		{
			p_alpha = alpha + j;
			p_edgePot = edgePot + maxState * (y[j+1] + maxState * j);
			sumPot = 0;
			for (int k = 0; k < nStates[j]; k++)
			{
				prob[k] = p_alpha[0] * p_edgePot[k];
				sumPot += prob[k];
				p_alpha += nNodes;
			}
			for (int k = 0; k < nStates[j]; k++)
				prob[k] /= sumPot;
			y[j] = sample(nStates[j], prob);
		}

		for (int j = 0; j < nNodes; j++)
			samples[i + size * j] = y[j] + 1;
	}

	UNPROTECT(7);
	return(_samples);
}
