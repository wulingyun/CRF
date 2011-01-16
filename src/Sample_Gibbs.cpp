#include "CRF.h"

SEXP Sample_Gibbs(SEXP _crf, SEXP _size, SEXP _burnIn, SEXP _start)
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

	SEXP _nAdj, _adjEdges, _temp;
	PROTECT(_nAdj = AS_INTEGER(getListElement(_crf, "n.adj")));
	PROTECT(_adjEdges = AS_LIST(getListElement(_crf, "adj.edges")));
	int *nAdj = INTEGER_POINTER(_nAdj);
	int **adjEdges = (int **) R_alloc(nNodes, sizeof(int *));
	for (int i = 0; i < nNodes; i++)
	{
		PROTECT(_temp = AS_INTEGER(VECTOR_ELT(_adjEdges, i)));
		adjEdges[i] = INTEGER_POINTER(_temp);
	}

	SEXP _nodePot, _edgePot;
	PROTECT(_nodePot = AS_NUMERIC(getListElement(_crf, "node.pot")));
	PROTECT(_edgePot = AS_NUMERIC(getListElement(_crf, "edge.pot")));
	double *nodePot = NUMERIC_POINTER(_nodePot);
	double *edgePot = NUMERIC_POINTER(_edgePot);

	PROTECT(_size = AS_INTEGER(_size));
	int size = INTEGER_POINTER(_size)[0];
	PROTECT(_burnIn = AS_INTEGER(_burnIn));
	int burnIn = INTEGER_POINTER(_burnIn)[0];
	PROTECT(_start = AS_INTEGER(_start));
	int *start = INTEGER_POINTER(_start);

	SEXP _samples;
	PROTECT(_samples = NEW_INTEGER(size * nNodes));
	setDim2(_samples, size, nNodes);
	int *samples = INTEGER_POINTER(_samples);
	setValues(_samples, samples, 0);

	int *y = (int *) R_alloc(nNodes, sizeof(int));
	for (int i = 0; i < nNodes; i++)
		y[i] = start[i] - 1;

	double sumProb, *prob = (double *) R_alloc(maxState, sizeof(double));

	int e, n, n1, n2;
	int *p_samples;
	double *p_nodePot, *p_edgePot;

	GetRNGstate();
	for (int iter = 0; iter < burnIn+size; iter++)
	{
		for (int i = 0; i < nNodes; i++)
		{
			n = nStates[i];
			p_nodePot = nodePot + i;
			for (int j = 0; j < n; j++)
			{
				prob[j] = p_nodePot[0];
				p_nodePot += nNodes;
			}
			for (int j = 0; j < nAdj[i]; j++)
			{
				e = adjEdges[i][j] - 1;
				n1 = edges[e] - 1;
				n2 = edges[e + nEdges] - 1;
				if (n1 == i)
				{
					p_edgePot = edgePot + maxState * (y[n2] + maxState * e);
					for (int k = 0; k < n; k++)
						prob[k] *= p_edgePot[k];
				}
				else
				{
					p_edgePot = edgePot + y[n1] + maxState * maxState * e;
					for (int k = 0; k < n; k++)
					{
						prob[k] *= p_edgePot[0];
						p_edgePot += maxState;
					}
				}
			}
			sumProb = 0;
			for (int j = 0; j < n; j++)
				sumProb += prob[j];
			for (int j = 0; j < n; j++)
				prob[j] /= sumProb;

			y[i] = sample(n, prob);
		}

		if (iter >= burnIn)
		{
			p_samples = samples + iter - burnIn;
			for (int i = 0; i < nNodes; i++)
			{
				p_samples[0] = y[i] + 1;
				p_samples += size;
			}
		}
	}
	PutRNGstate();

	UNPROTECT(13 + nNodes);

	return(_samples);
}
