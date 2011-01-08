#include <time.h>
#include "CRF.h"

SEXP Sample_Exact(SEXP _crf, SEXP _size)
{
	SEXP _nNodes, _nEdges, _edges, _nStates, _maxState;
	PROTECT(_nNodes = AS_INTEGER(getListElement(_crf, "n.nodes")));
	PROTECT(_nEdges = AS_INTEGER(getListElement(_crf, "n.edges")));
	PROTECT(_edges = AS_INTEGER(getListElement(_crf, "edges")));
	PROTECT(_nStates = AS_INTEGER(getListElement(_crf, "n.states")));
	PROTECT(_maxState = AS_INTEGER(getListElement(_crf, "max.state")));
	int nNodes = INTEGER_POINTER(_nNodes)[0];
	int nEdges = INTEGER_POINTER(_nEdges)[0];
	int *edges = INTEGER_POINTER(_edges);
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

	double pot, Z = 0;
	int index;
	while(1)
	{
		pot = 1;

		/* Node potentials */
		for (int i = 0; i < nNodes; i++)
			pot *= nodePot[i + nNodes * y[i]];

		/* Edge potentials */
		for (int i = 0; i < nEdges; i++)
			pot *= edgePot[y[edges[i]-1] + maxState * (y[edges[i+nEdges]-1] + maxState * i)];

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

	/* Sampling */

	double cutoff, cumulativePot;
	GetRNGstate();
	for (int k = 0; k < size; k++)
	{
		for (int i = 0; i < nNodes; i++)
			y[i] = 0;

		cutoff = unif_rand();
		cumulativePot = 0;
		while(1)
		{
			pot = 1;

			/* Node potentials */
			for (int i = 0; i < nNodes; i++)
				pot *= nodePot[i + nNodes * y[i]];

			/* Edge potentials */
			for (int i = 0; i < nEdges; i++)
				pot *= edgePot[y[edges[i]-1] + maxState * (y[edges[i+nEdges]-1] + maxState * i)];

			/* Update cumulative potential */
			cumulativePot += pot;

			if (cumulativePot/Z > cutoff)
				break;

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

		for (int i = 0; i < nNodes; i++)
			samples[k + size * i] = y[i] + 1;
	}
	PutRNGstate();

	UNPROTECT(9);
	return(_samples);
}
