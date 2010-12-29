#include "CRF.h"

SEXP Decode_Chain(SEXP _crf)
{
	int nNodes, nEdges, *edges, *nStates, maxState, *labels;
	SEXP _nNodes, _nEdges, _edges, _nStates, _maxState, _labels;
	double *nodePot, *edgePot;
	SEXP _nodePot, _edgePot;

	PROTECT(_nNodes = AS_INTEGER(getListElement(_crf, "n.nodes")));
	PROTECT(_nEdges = AS_INTEGER(getListElement(_crf, "n.edges")));
	PROTECT(_edges = AS_INTEGER(getListElement(_crf, "edges")));
	PROTECT(_nStates = AS_INTEGER(getListElement(_crf, "n.states")));
	PROTECT(_maxState = AS_INTEGER(getListElement(_crf, "max.state")));
	nNodes = INTEGER_POINTER(_nNodes)[0];
	nEdges = INTEGER_POINTER(_nEdges)[0];
	edges = INTEGER_POINTER(_edges);
	nStates = INTEGER_POINTER(_nStates);
	maxState = INTEGER_POINTER(_maxState)[0];

	PROTECT(_nodePot = AS_NUMERIC(getListElement(_crf, "node.pot")));
	PROTECT(_edgePot = AS_NUMERIC(getListElement(_crf, "edge.pot")));
	nodePot = NUMERIC_POINTER(_nodePot);
	edgePot = NUMERIC_POINTER(_edgePot);

	PROTECT(_labels = NEW_INTEGER(nNodes));
	labels = INTEGER_POINTER(_labels);

	int *y = (int *) R_alloc(nNodes, sizeof(int));
	for (int i = 0; i < nNodes; i++)
		y[i] = 0;

	double pot, maxPot = -1;
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

		/* Record the best potentials */
		if (pot > maxPot)
		{
			maxPot = pot;
			for (int i = 0; i < nNodes; i++)
				labels[i] = y[i] + 1;
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

		if (index == nNodes)
			break;
	}

	UNPROTECT(8);
	return(_labels);
}
