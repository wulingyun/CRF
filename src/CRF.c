#include "CRF.h"

void openCRF(CRFinfo *crf, SEXP _crf)
{
	PROTECT(crf->_nNodes = AS_INTEGER(getListElement(_crf, "n.nodes")));
	PROTECT(crf->_nEdges = AS_INTEGER(getListElement(_crf, "n.edges")));
	PROTECT(crf->_edges = AS_INTEGER(getListElement(_crf, "edges")));
	PROTECT(crf->_nStates = AS_INTEGER(getListElement(_crf, "n.states")));
	PROTECT(crf->_maxState = AS_INTEGER(getListElement(_crf, "max.state")));
	crf->nNodes = INTEGER_POINTER(crf->_nNodes)[0];
	crf->nEdges = INTEGER_POINTER(crf->_nEdges)[0];
	crf->edges = INTEGER_POINTER(crf->_edges);
	crf->nStates = INTEGER_POINTER(crf->_nStates);
	crf->maxState = INTEGER_POINTER(crf->_maxState)[0];

	PROTECT(crf->_nAdj = AS_INTEGER(getListElement(_crf, "n.adj")));
	PROTECT(crf->_adjNodes = AS_LIST(getListElement(_crf, "adj.nodes")));
	PROTECT(crf->_adjEdges = AS_LIST(getListElement(_crf, "adj.edges")));
	crf->nAdj = INTEGER_POINTER(crf->_nAdj);
	crf->adjNodes = (int **) R_alloc(crf->nNodes, sizeof(int *));
	crf->adjEdges = (int **) R_alloc(crf->nNodes, sizeof(int *));
	SEXP _temp;
	for (int i = 0; i < crf->nNodes; i++)
	{
		PROTECT(_temp = AS_INTEGER(VECTOR_ELT(crf->_adjNodes, i)));
		crf->adjNodes[i] = INTEGER_POINTER(_temp);
		PROTECT(_temp = AS_INTEGER(VECTOR_ELT(crf->_adjEdges, i)));
		crf->adjEdges[i] = INTEGER_POINTER(_temp);
	}

	PROTECT(crf->_nodePot = AS_NUMERIC(getListElement(_crf, "node.pot")));
	PROTECT(crf->_edgePot = AS_NUMERIC(getListElement(_crf, "edge.pot")));
	crf->nodePot = NUMERIC_POINTER(crf->_nodePot);
	crf->edgePot = NUMERIC_POINTER(crf->_edgePot);
}

void closeCRF(CRFinfo *crf)
{
	UNPROTECT(10 + crf->nNodes * 2);
}
