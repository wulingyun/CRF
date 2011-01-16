#include "CRF.h"

CRF::CRF(SEXP _crf)
{
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

	PROTECT(_nAdj = AS_INTEGER(getListElement(_crf, "n.adj")));
	PROTECT(_adjNodes = AS_LIST(getListElement(_crf, "adj.nodes")));
	PROTECT(_adjEdges = AS_LIST(getListElement(_crf, "adj.edges")));
	nAdj = INTEGER_POINTER(_nAdj);
	adjNodes = (int **) R_alloc(nNodes, sizeof(int *));
	adjEdges = (int **) R_alloc(nNodes, sizeof(int *));
	SEXP _temp;
	for (int i = 0; i < nNodes; i++)
	{
		PROTECT(_temp = AS_INTEGER(VECTOR_ELT(_adjNodes, i)));
		adjNodes[i] = INTEGER_POINTER(_temp);
		PROTECT(_temp = AS_INTEGER(VECTOR_ELT(_adjEdges, i)));
		adjEdges[i] = INTEGER_POINTER(_temp);
	}

	PROTECT(_nodePot = AS_NUMERIC(getListElement(_crf, "node.pot")));
	PROTECT(_edgePot = AS_NUMERIC(getListElement(_crf, "edge.pot")));
	nodePot = NUMERIC_POINTER(_nodePot);
	edgePot = NUMERIC_POINTER(_edgePot);
}

CRF::~CRF()
{
	UNPROTECT(10 + nNodes * 2);
}
