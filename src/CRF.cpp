#include "CRF.h"

const R_CallMethodDef callMethods[] = {
  {"Decode_Exact", (DL_FUNC) &Decode_Exact, 1},
  {"Decode_Chain", (DL_FUNC) &Decode_Chain, 1},
  {"Decode_Tree", (DL_FUNC) &Decode_Tree, 1},
  {"Decode_Cutset", (DL_FUNC) &Decode_Cutset, 3},
  {"Decode_Junction", (DL_FUNC) &Decode_Junction, 1},
  {"Decode_Sample", (DL_FUNC) &Decode_Sample, 2},
  {"Decode_LBP", (DL_FUNC) &Decode_LBP, 4},
  {"Decode_TRBP", (DL_FUNC) &Decode_TRBP, 4},
  {"Decode_Greedy", (DL_FUNC) &Decode_Greedy, 3},
  {"Decode_ICM", (DL_FUNC) &Decode_ICM, 3},
  {"Infer_Exact", (DL_FUNC) &Infer_Exact, 1},
  {"Infer_Chain", (DL_FUNC) &Infer_Chain, 1},
  {"Infer_Tree", (DL_FUNC) &Infer_Tree, 1},
  {"Infer_Cutset", (DL_FUNC) &Infer_Cutset, 2},
  {"Infer_Junction", (DL_FUNC) &Infer_Junction, 1},
  {"Infer_Sample", (DL_FUNC) &Infer_Sample, 2},
  {"Infer_LBP", (DL_FUNC) &Infer_LBP, 4},
  {"Infer_TRBP", (DL_FUNC) &Infer_TRBP, 4},
  {"Sample_Exact", (DL_FUNC) &Sample_Exact, 2},
  {"Sample_Chain", (DL_FUNC) &Sample_Chain, 2},
  {"Sample_Tree", (DL_FUNC) &Sample_Tree, 2},
  {"Sample_Cutset", (DL_FUNC) &Sample_Cutset, 3},
  {"Sample_Junction", (DL_FUNC) &Sample_Junction, 2},
  {"Sample_Gibbs", (DL_FUNC) &Sample_Gibbs, 4},
  {"Make_AdjInfo", (DL_FUNC) &Make_AdjInfo, 1},
  {"MRF_Update", (DL_FUNC) &MRF_Update, 1},
  {"CRF_Update", (DL_FUNC) &CRF_Update, 5},
  {"MRF_Stat", (DL_FUNC) &MRF_Stat, 2},
  {"MRF_NLL", (DL_FUNC) &MRF_NLL, 5},
  {"CRF_NLL", (DL_FUNC) &CRF_NLL, 9},
  {"Clamp_Reset", (DL_FUNC) &Clamp_Reset, 1},
  {"Get_Potential", (DL_FUNC) &Get_Potential, 2},
  {"Get_LogPotential", (DL_FUNC) &Get_LogPotential, 2},
  {"Calc_Frequency", (DL_FUNC) &Calc_Frequency, 2},
  {NULL, NULL, 0}
};

void R_init_CRF(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

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

void CRF::Set_Data(SEXP _newcrf)
{
	_crf = _newcrf;
	PROTECT(_nNodes = AS_INTEGER(GetVar(_crf, "n.nodes")));
	PROTECT(_nEdges = AS_INTEGER(GetVar(_crf, "n.edges")));
	PROTECT(_edges = AS_INTEGER(GetVar(_crf, "edges")));
	PROTECT(_nStates = AS_INTEGER(GetVar(_crf, "n.states")));
	PROTECT(_maxState = AS_INTEGER(GetVar(_crf, "max.state")));
	nNodes = INTEGER_POINTER(_nNodes)[0];
	nEdges = INTEGER_POINTER(_nEdges)[0];
	edges = INTEGER_POINTER(_edges);
	nStates = INTEGER_POINTER(_nStates);
	maxState = INTEGER_POINTER(_maxState)[0];

	PROTECT(_nAdj = AS_INTEGER(GetVar(_crf, "n.adj")));
	PROTECT(_adjNodes = AS_LIST(GetVar(_crf, "adj.nodes")));
	PROTECT(_adjEdges = AS_LIST(GetVar(_crf, "adj.edges")));
	nAdj = INTEGER_POINTER(_nAdj);
	adjNodes = (int **) R_alloc(nNodes, sizeof(int *));
	adjEdges = (int **) R_alloc(nNodes, sizeof(int *));
	SEXP _temp;
	for (int i = 0; i < nNodes; i++)
	{
	  SET_VECTOR_ELT(_adjNodes, i, _temp = AS_INTEGER(VECTOR_ELT(_adjNodes, i)));
		adjNodes[i] = INTEGER_POINTER(_temp);
		SET_VECTOR_ELT(_adjEdges, i, _temp = AS_INTEGER(VECTOR_ELT(_adjEdges, i)));
		adjEdges[i] = INTEGER_POINTER(_temp);
	}

	PROTECT(_nodePot = AS_NUMERIC(GetVar(_crf, "node.pot")));
	PROTECT(_edgePot = AS_LIST(GetVar(_crf, "edge.pot")));
	nodePot = NUMERIC_POINTER(_nodePot);
	edgePot = (double **) R_alloc(nEdges, sizeof(double *));
	nEdgeStates = (int *) R_alloc(nEdges, sizeof(int));
	for (int i = 0; i < nEdges; i++)
	{
	  SET_VECTOR_ELT(_edgePot, i, _temp = AS_NUMERIC(VECTOR_ELT(_edgePot, i)));
		edgePot[i] = NUMERIC_POINTER(_temp);
		nEdgeStates[i] = nStates[EdgesBegin(i)] * nStates[EdgesEnd(i)];
	}

	numProtect = 10;
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
