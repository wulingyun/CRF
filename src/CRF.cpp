#include "CRF.h"

CRF::CRF()
{
	nNodes = 0;
	nEdges = 0;
	edges = NULL;
	nStates = NULL;
	maxState = 0;

	nAdj = NULL;
	adjNodes = NULL;
	adjEdges = NULL;

	nodePot = NULL;
	edgePot = NULL;
	nEdgeStates = NULL;

	labels = NULL;

	nodeBel = NULL;
	edgeBel = NULL;
	logZ = NULL;

	samples = NULL;
	nSamples = 0;

	maxNodePot = NULL;
	maxEdgePot = NULL;
	messages = NULL;

	numProtect = 0;
}

CRF::CRF(SEXP _crf)
{
	Set_Data(_crf);
}

CRF::~CRF()
{
	UNPROTECT(numProtect);
}

void CRF::Set_Data(SEXP _crf)
{
	PROTECT(_nNodes = AS_INTEGER(GetListElement(_crf, "n.nodes")));
	PROTECT(_nEdges = AS_INTEGER(GetListElement(_crf, "n.edges")));
	PROTECT(_edges = AS_INTEGER(GetListElement(_crf, "edges")));
	PROTECT(_nStates = AS_INTEGER(GetListElement(_crf, "n.states")));
	PROTECT(_maxState = AS_INTEGER(GetListElement(_crf, "max.state")));
	nNodes = INTEGER_POINTER(_nNodes)[0];
	nEdges = INTEGER_POINTER(_nEdges)[0];
	edges = INTEGER_POINTER(_edges);
	nStates = INTEGER_POINTER(_nStates);
	maxState = INTEGER_POINTER(_maxState)[0];

	PROTECT(_nAdj = AS_INTEGER(GetListElement(_crf, "n.adj")));
	PROTECT(_adjNodes = AS_LIST(GetListElement(_crf, "adj.nodes")));
	PROTECT(_adjEdges = AS_LIST(GetListElement(_crf, "adj.edges")));
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

	PROTECT(_nodePot = AS_NUMERIC(GetListElement(_crf, "node.pot")));
	PROTECT(_edgePot = AS_LIST(GetListElement(_crf, "edge.pot")));
	nodePot = NUMERIC_POINTER(_nodePot);
	edgePot = (double **) R_alloc(nEdges, sizeof(double *));
	nEdgeStates = (int *) R_alloc(nEdges, sizeof(int));
	for (int i = 0; i < nEdges; i++)
	{
		PROTECT(_temp = AS_NUMERIC(VECTOR_ELT(_edgePot, i)));
		edgePot[i] = NUMERIC_POINTER(_temp);
		nEdgeStates[i] = nStates[EdgesBegin(i)] * nStates[EdgesEnd(i)];
	}

	numProtect = 10 + nNodes * 2 + nEdges;
}

void CRF::Init_Labels()
{
	PROTECT(_labels = NEW_INTEGER(nNodes));
	labels = INTEGER_POINTER(_labels);
	SetValues(_labels, labels, 1);
	numProtect++;
}

void CRF::Init_NodeBel()
{
	PROTECT(_nodeBel = NEW_NUMERIC(nNodes * maxState));
	SetDim2(_nodeBel, nNodes, maxState);
	nodeBel = NUMERIC_POINTER(_nodeBel);
	SetValues(_nodeBel, nodeBel, 0.0);
	numProtect++;
}

void CRF::Init_EdgeBel()
{
	PROTECT(_edgeBel = NEW_LIST(nEdges));
	edgeBel = (double **) R_alloc(nEdges, sizeof(double *));
	SEXP _temp;
	for (int i = 0; i < nEdges; i++)
	{
		PROTECT(_temp = NEW_NUMERIC(nEdgeStates[i]));
		SetDim2(_temp, nStates[EdgesBegin(i)], nStates[EdgesEnd(i)]);
		edgeBel[i] = NUMERIC_POINTER(_temp);
		SetValues(_temp, edgeBel[i], 0.0);
		SET_VECTOR_ELT(_edgeBel, i, _temp);
	}
	numProtect += nEdges + 1;
}

void CRF::Init_LogZ()
{
	PROTECT(_logZ = NEW_NUMERIC(1));
	logZ = NUMERIC_POINTER(_logZ);
	*logZ = 0;
	numProtect++;
}

void CRF::Init_Belief()
{
	Init_NodeBel();
	Init_EdgeBel();
	Init_LogZ();
	PROTECT(_belief = NEW_LIST(3));
	SetListElement(_belief, 0, "node.bel", _nodeBel);
	SetListElement(_belief, 1, "edge.bel", _edgeBel);
	SetListElement(_belief, 2, "logZ", _logZ);
	numProtect++;
}

void CRF::Init_Samples(int size)
{
	nSamples = size;
	PROTECT(_samples = NEW_INTEGER(size * nNodes));
	SetDim2(_samples, size, nNodes);
	samples = INTEGER_POINTER(_samples);
	SetValues(_samples, samples, 0);
	numProtect++;
}

void CRF::Init_Samples(SEXP _size)
{
	Init_Samples(INTEGER_POINTER(AS_INTEGER(_size))[0]);
}

void CRF::Set_Samples(SEXP _otherSamples)
{
	PROTECT(_samples = AS_INTEGER(_otherSamples));
	samples = INTEGER_POINTER(_samples);
	nSamples = length(_samples) / nNodes;
	numProtect++;
}

void CRF::Normalize_NodeBel()
{
	double sumBel;
	for (int i = 0; i < nNodes; i++)
	{
		sumBel = 0;
		for (int j = 0; j < nStates[i]; j++)
			sumBel += NodeBel(i, j);
		for (int j = 0; j < nStates[i]; j++)
			NodeBel(i, j) /= sumBel;
	}
}

void CRF::Normalize_EdgeBel()
{
	int n1, n2;
	double sumBel;
	for (int i = 0; i < nEdges; i++)
	{
		n1 = EdgesBegin(i);
		n2 = EdgesEnd(i);
		sumBel = 0;
		for (int j = 0; j < nStates[n2]; j++)
		{
			for (int k = 0; k < nStates[n1]; k++)
				sumBel += EdgeBel(i, k, j);
		}
		for (int j = 0; j < nStates[n2]; j++)
		{
			for (int k = 0; k < nStates[n1]; k++)
				EdgeBel(i, k, j) /= sumBel;
		}
	}
}
