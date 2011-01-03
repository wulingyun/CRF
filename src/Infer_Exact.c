#include "CRF.h"

SEXP Infer_Exact(SEXP _crf)
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

	SEXP _nodeBel, _edgeBel, _logZ;
	PROTECT(_nodeBel = NEW_NUMERIC(nNodes * maxState));
	PROTECT(_edgeBel = NEW_NUMERIC(maxState * maxState * nEdges));
	PROTECT(_logZ = NEW_NUMERIC(1));
	setDim2(_nodeBel, nNodes, maxState);
	setDim3(_edgeBel, maxState, maxState, nEdges);
	double *nodeBel = NUMERIC_POINTER(_nodeBel);
	double *edgeBel = NUMERIC_POINTER(_edgeBel);
	double *logZ = NUMERIC_POINTER(_logZ);

	for (int i = 0; i < length(_nodeBel); i++)
		nodeBel[i] = 0;
	for (int i = 0; i < length(_edgeBel); i++)
		edgeBel[i] = 0;
	*logZ = 0;

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

		/* Node belief */
		for (int i = 0; i < nNodes; i++)
			nodeBel[i + nNodes * y[i]] += pot;

		/* Edge belief */
		for (int i = 0; i < nEdges; i++)
			edgeBel[y[edges[i]-1] + maxState * (y[edges[i+nEdges]-1] + maxState * i)] += pot;

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

	/* Normalization */
	for (int i = 0; i < length(_nodeBel); i++)
		nodeBel[i] /= Z;
	for (int i = 0; i < length(_edgeBel); i++)
		edgeBel[i] /= Z;
	*logZ = log(Z);

	SEXP _belief;
	PROTECT(_belief = NEW_LIST(3));
	setListElement(_belief, 0, "node.bel", _nodeBel);
	setListElement(_belief, 1, "edge.bel", _edgeBel);
	setListElement(_belief, 2, "logZ", _logZ);

	UNPROTECT(11);
	return(_belief);
}
