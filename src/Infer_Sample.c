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

	SEXP _nodePot, _edgePot;
	PROTECT(_nodePot = AS_NUMERIC(getListElement(_crf, "node.pot")));
	PROTECT(_edgePot = AS_NUMERIC(getListElement(_crf, "edge.pot")));
	double *nodePot = NUMERIC_POINTER(_nodePot);
	double *edgePot = NUMERIC_POINTER(_edgePot);

	PROTECT(_samples = AS_INTEGER(_samples));
	int *samples = INTEGER_POINTER(_samples);
	int nSamples = length(_samples) / nNodes;

	SEXP _nodeBel, _edgeBel, _logZ;
	PROTECT(_nodeBel = NEW_NUMERIC(nNodes * maxState));
	PROTECT(_edgeBel = NEW_NUMERIC(maxState * maxState * nEdges));
	PROTECT(_logZ = NEW_NUMERIC(1));
	setDim2(_nodeBel, nNodes, maxState);
	setDim3(_edgeBel, maxState, maxState, nEdges);
	double *nodeBel = NUMERIC_POINTER(_nodeBel);
	double *edgeBel = NUMERIC_POINTER(_edgeBel);
	double *logZ = NUMERIC_POINTER(_logZ);
	setValues(_nodeBel, nodeBel, 0);
	setValues(_edgeBel, edgeBel, 0);
	*logZ = 0;

	double pot, maxPot = -1;
	int k, maxSample = -1;
	for (int i = 0; i < nSamples; i++)
	{
		pot = 1;
		for (int j = 0; j < nNodes; j++)
		{
			k = j + nNodes * (samples[i + nSamples * j] - 1);
			nodeBel[k]++;
			pot *= nodePot[k];
		}
		for (int j = 0; j < nEdges; j++)
		{
			k = samples[i + nSamples * (edges[j] - 1)] - 1 + maxState * (samples[i + nSamples * (edges[j + nEdges] - 1)] - 1 + maxState * j);
			edgeBel[k]++;
			pot *= edgePot[k];
		}

		if (pot > maxPot)
		{
			maxPot = pot;
			maxSample = i;
		}
	}
	int same, freq = 0;
	for (int i = 0; i < nSamples; i++)
	{
		same = 1;
		for (int j = 0; j < nNodes; j++)
		{
			if (samples[i + nSamples * j] != samples[maxSample + nSamples * j])
			{
				same = 0;
				break;
			}
		}
		if (same)
			freq++;
	}
	for (int i = 0; i < length(_nodeBel); i++)
		nodeBel[i] /= nSamples;
	for (int i = 0; i < length(_edgeBel); i++)
		edgeBel[i] /= nSamples;
	*logZ = log(maxPot * nSamples / freq);

	SEXP _belief;
	PROTECT(_belief = NEW_LIST(3));
	setListElement(_belief, 0, "node.bel", _nodeBel);
	setListElement(_belief, 1, "edge.bel", _edgeBel);
	setListElement(_belief, 2, "logZ", _logZ);

	UNPROTECT(11);
	return(_belief);
}
