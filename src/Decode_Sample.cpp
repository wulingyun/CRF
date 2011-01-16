#include "CRF.h"

SEXP Decode_Sample(SEXP _crf, SEXP _samples)
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

	PROTECT(_samples = AS_INTEGER(_samples));
	int *samples = INTEGER_POINTER(_samples);
	int nSamples = length(_samples) / nNodes;

	SEXP _labels;
	PROTECT(_labels = NEW_INTEGER(nNodes));
	int *labels = INTEGER_POINTER(_labels);
	setValues(_labels, labels, -1);

	double pot, maxPot = -1;
	int k, maxSample = -1;
	for (int i = 0; i < nSamples; i++)
	{
		pot = 1;
		for (int j = 0; j < nNodes; j++)
		{
			k = j + nNodes * (samples[i + nSamples * j] - 1);
			pot *= nodePot[k];
		}
		for (int j = 0; j < nEdges; j++)
		{
			k = samples[i + nSamples * (edges[j] - 1)] - 1 + maxState * (samples[i + nSamples * (edges[j + nEdges] - 1)] - 1 + maxState * j);
			pot *= edgePot[k];
		}

		if (pot > maxPot)
		{
			maxPot = pot;
			maxSample = i;
		}
	}

	for (int i = 0; i < nNodes; i++)
		labels[i] = samples[maxSample + nSamples * i];

	UNPROTECT(8);

	return(_labels);
}
