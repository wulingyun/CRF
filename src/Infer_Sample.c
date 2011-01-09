#include "CRF.h"

SEXP Infer_Sample(SEXP _crf, SEXP _samples)
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

	PROTECT(_samples = AS_INTEGER(_samples));
	int *samples = INTEGER_POINTER(_samples);
	int nSamples = length(_samples) / nNodes;

	SEXP _nodeBel, _edgeBel;
	PROTECT(_nodeBel = NEW_NUMERIC(nNodes * maxState));
	PROTECT(_edgeBel = NEW_NUMERIC(maxState * maxState * nEdges));
	setDim2(_nodeBel, nNodes, maxState);
	setDim3(_edgeBel, maxState, maxState, nEdges);
	double *nodeBel = NUMERIC_POINTER(_nodeBel);
	double *edgeBel = NUMERIC_POINTER(_edgeBel);
	setValues(_nodeBel, nodeBel, 0);
	setValues(_edgeBel, edgeBel, 0);

	for (int i = 0; i < nSamples; i++)
	{
		for (int j = 0; j < nNodes; j++)
			nodeBel[j + nNodes * (samples[i + nSamples * j] - 1)]++;
		for (int j = 0; j < nEdges; j++)
			edgeBel[samples[i + nSamples * (edges[j] - 1)] - 1 + maxState * (samples[i + nSamples * (edges[j + nEdges] - 1)] - 1 + maxState * j)]++;
	}
	for (int i = 0; i < length(_nodeBel); i++)
		nodeBel[i] /= nSamples;
	for (int i = 0; i < length(_edgeBel); i++)
		edgeBel[i] /= nSamples;

	SEXP _belief;
	PROTECT(_belief = NEW_LIST(2));
	setListElement(_belief, 0, "node.bel", _nodeBel);
	setListElement(_belief, 1, "edge.bel", _edgeBel);

	UNPROTECT(8);
	return(_belief);
}
