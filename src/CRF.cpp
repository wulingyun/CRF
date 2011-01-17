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

	numProtect = 10 + nNodes * 2;
}

CRF::~CRF()
{
	UNPROTECT(numProtect);
}

void CRF::Init_Decoding()
{
	PROTECT(_labels = NEW_INTEGER(nNodes));
	labels = INTEGER_POINTER(_labels);
	setValues(_labels, labels, -1);
	numProtect++;
}

void CRF::Init_Decoding2()
{
	Init_Decoding();
	PROTECT(_nodeBel = NEW_NUMERIC(nNodes * maxState));
	setDim2(_nodeBel, nNodes, maxState);
	nodeBel = NUMERIC_POINTER(_nodeBel);
	setValues(_nodeBel, nodeBel, 0);
	numProtect++;
}

void CRF::Init_Inference()
{
	PROTECT(_nodeBel = NEW_NUMERIC(nNodes * maxState));
	PROTECT(_edgeBel = NEW_NUMERIC(maxState * maxState * nEdges));
	PROTECT(_logZ = NEW_NUMERIC(1));
	PROTECT(_belief = NEW_LIST(3));
	setDim2(_nodeBel, nNodes, maxState);
	setDim3(_edgeBel, maxState, maxState, nEdges);
	setListElement(_belief, 0, "node.bel", _nodeBel);
	setListElement(_belief, 1, "edge.bel", _edgeBel);
	setListElement(_belief, 2, "logZ", _logZ);
	nodeBel = NUMERIC_POINTER(_nodeBel);
	edgeBel = NUMERIC_POINTER(_edgeBel);
	logZ = NUMERIC_POINTER(_logZ);
	setValues(_nodeBel, nodeBel, 0);
	setValues(_edgeBel, edgeBel, 0);
	*logZ = 0;
	numProtect += 4;
}

void CRF::Init_Sampling(int size)
{
	PROTECT(_samples = NEW_INTEGER(size * nNodes));
	setDim2(_samples, size, nNodes);
	samples = INTEGER_POINTER(_samples);
	setValues(_samples, samples, 0);
	numProtect++;
}

void CRF::Init_Sampling2(int size)
{
	Init_Sampling(size);
	PROTECT(_nodeBel = NEW_NUMERIC(nNodes * maxState));
	PROTECT(_edgeBel = NEW_NUMERIC(maxState * maxState * nEdges));
	setDim2(_nodeBel, nNodes, maxState);
	setDim3(_edgeBel, maxState, maxState, nEdges);
	nodeBel = NUMERIC_POINTER(_nodeBel);
	edgeBel = NUMERIC_POINTER(_edgeBel);
	setValues(_nodeBel, nodeBel, 0);
	setValues(_edgeBel, edgeBel, 0);
	numProtect += 2;
}
