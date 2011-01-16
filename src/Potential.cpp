#include "CRF.h"

SEXP Potential(SEXP _crf, SEXP _configuration)
{
	SEXP _nNodes, _nEdges, _edges, _maxState;
	PROTECT(_nNodes = AS_INTEGER(getListElement(_crf, "n.nodes")));
	PROTECT(_nEdges = AS_INTEGER(getListElement(_crf, "n.edges")));
	PROTECT(_edges = AS_INTEGER(getListElement(_crf, "edges")));
	PROTECT(_maxState = AS_INTEGER(getListElement(_crf, "max.state")));
	int nNodes = INTEGER_POINTER(_nNodes)[0];
	int nEdges = INTEGER_POINTER(_nEdges)[0];
	int *edges = INTEGER_POINTER(_edges);
	int maxState = INTEGER_POINTER(_maxState)[0];

	SEXP _nodePot, _edgePot;
	PROTECT(_nodePot = AS_NUMERIC(getListElement(_crf, "node.pot")));
	PROTECT(_edgePot = AS_NUMERIC(getListElement(_crf, "edge.pot")));
	double *nodePot = NUMERIC_POINTER(_nodePot);
	double *edgePot = NUMERIC_POINTER(_edgePot);

	PROTECT(_configuration = AS_INTEGER(_configuration));
	int *configuration = INTEGER_POINTER(_configuration);

	SEXP _potential;
	PROTECT(_potential = NEW_NUMERIC(1));
	double *potential = NUMERIC_POINTER(_potential);

	int *y = (int *) R_alloc(nNodes, sizeof(int));
	for (int i = 0; i < nNodes; i++)
		y[i] = configuration[i] - 1;

	*potential = 1;

	/* Node potentials */
	for (int i = 0; i < nNodes; i++)
		*potential *= nodePot[i + nNodes * y[i]];

	/* Edge potentials */
	for (int i = 0; i < nEdges; i++)
		*potential *= edgePot[y[edges[i]-1] + maxState * (y[edges[i+nEdges]-1] + maxState * i)];

	UNPROTECT(8);
	return(_potential);
}

SEXP Log_Potential(SEXP _crf, SEXP _configuration)
{
	SEXP _nNodes, _nEdges, _edges, _maxState;
	PROTECT(_nNodes = AS_INTEGER(getListElement(_crf, "n.nodes")));
	PROTECT(_nEdges = AS_INTEGER(getListElement(_crf, "n.edges")));
	PROTECT(_edges = AS_INTEGER(getListElement(_crf, "edges")));
	PROTECT(_maxState = AS_INTEGER(getListElement(_crf, "max.state")));
	int nNodes = INTEGER_POINTER(_nNodes)[0];
	int nEdges = INTEGER_POINTER(_nEdges)[0];
	int *edges = INTEGER_POINTER(_edges);
	int maxState = INTEGER_POINTER(_maxState)[0];

	SEXP _nodePot, _edgePot;
	PROTECT(_nodePot = AS_NUMERIC(getListElement(_crf, "node.pot")));
	PROTECT(_edgePot = AS_NUMERIC(getListElement(_crf, "edge.pot")));
	double *nodePot = NUMERIC_POINTER(_nodePot);
	double *edgePot = NUMERIC_POINTER(_edgePot);

	PROTECT(_configuration = AS_INTEGER(_configuration));
	int *configuration = INTEGER_POINTER(_configuration);

	SEXP _potential;
	PROTECT(_potential = NEW_NUMERIC(1));
	double *potential = NUMERIC_POINTER(_potential);

	int *y = (int *) R_alloc(nNodes, sizeof(int));
	for (int i = 0; i < nNodes; i++)
		y[i] = configuration[i] - 1;

	*potential = 0;

	/* Node potentials */
	for (int i = 0; i < nNodes; i++)
		*potential += log(nodePot[i + nNodes * y[i]]);

	/* Edge potentials */
	for (int i = 0; i < nEdges; i++)
		*potential += log(edgePot[y[edges[i]-1] + maxState * (y[edges[i+nEdges]-1] + maxState * i)]);

	UNPROTECT(8);
	return(_potential);
}
