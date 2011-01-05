#include "CRF.h"

SEXP Decode_Tree(SEXP _crf)
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

	SEXP _nAdj, _adjNodes, _adjEdges;
	PROTECT(_nAdj = AS_INTEGER(getListElement(_crf, "n.adj")));
	PROTECT(_adjNodes = AS_INTEGER(getListElement(_crf, "adj.nodes")));
	PROTECT(_adjEdges = AS_INTEGER(getListElement(_crf, "adj.edges")));
	int *nAdj = INTEGER_POINTER(_nAdj);
	int **adjNodes = (int **) R_alloc(nNodes, sizeof(int *));
	int **adjEdges = (int **) R_alloc(nNodes, sizeof(int *));
	for (int i = 0; i < nNodes; i++)
	{
		adjNodes[i] = INTEGER_POINTER(VECTOR_ELT(_adjNodes, i));
		adjEdges[i] = INTEGER_POINTER(VECTOR_ELT(_adjEdges, i));
	}

	SEXP _nodePot, _edgePot;
	PROTECT(_nodePot = AS_NUMERIC(getListElement(_crf, "node.pot")));
	PROTECT(_edgePot = AS_NUMERIC(getListElement(_crf, "edge.pot")));
	double *nodePot = NUMERIC_POINTER(_nodePot);
	double *edgePot = NUMERIC_POINTER(_edgePot);

	SEXP _labels;
	PROTECT(_labels = NEW_INTEGER(nNodes));
	int *labels = INTEGER_POINTER(_labels);
	setValues(_labels, labels, -1);

	/* Tree BP */

	double *nodeBel = (double *) R_alloc(nNodes * maxState, sizeof(double));
	for (int i = 0; i < nNodes * maxState; i++)
		nodeBel[i] = nodePot[i];
	double *messages_1 = (double *) R_alloc(maxState * nEdges, sizeof(double));
	double *messages_2 = (double *) R_alloc(maxState * nEdges, sizeof(double));
	for (int i = 0; i < maxState * nEdges; i++)
		messages_1[i] = messages_2[i] = 0;
	int *nWaiting = (int *) R_alloc(nNodes, sizeof(int));
	int **waiting = (int **) R_alloc(nNodes, sizeof(int *));
	int *nUnsent = (int *) R_alloc(nNodes, sizeof(int));
	int **unsent = (int **) R_alloc(nNodes, sizeof(int *));
	for (int i = 0; i < nNodes; i++)
	{
		nWaiting[i] = nUnsent[i] = nAdj[i];
		waiting[i] = (int *) R_alloc(nAdj[i], sizeof(int));
		unsent[i] = (int *) R_alloc(nAdj[i], sizeof(int));
		for (int j = 0; j < nAdj[i]; j++)
			waiting[i][j] = unsent[i][j] = 1;
	}

	int nQueue;
	int *queue = (int *) R_alloc(nNodes, sizeof(int *));

	int n, r, e;
	int done = 0;
	while (!done)
	{
		done = 1;
		for (int i = 0; i < nNodes; i++)
		{
			if (nUnsent[i] == 0 || nWaiting[i] > 1)
				continue;

			nQueue = 0;
			if (nWaiting[i] == 1)
				for (int j = 0; j < nAdj[i]; j++)
					if (waiting[i][j] && unsent[i][j])
						queue[nQueue++] = j;
			else
				for (int j = 0; j < nAdj[i]; j++)
					if (unsent[i][j])
						queue[nQueue++] = j;

			if (nQueue > 0)
			{
				for (int j = 0; j < nQueue; j++)
				{
					n = queue[j];
					unsent[i][n] = 0;
					nUnsent[i]--;
					r = adjNodes[i][n];
					e = adjEdges[i][n];

					/* todo: calculate messages */

					for (int k = 0; k < nAdj[r]; k++)
						if (adjNodes[r][k] == i)
						{
							waiting[r][k] = 0;
							nWaiting[r]--;
						}
				}
				done = 0;
			}
		}
	}


	/* Get decode */

	double maxBel;
	double *p_nodeBel;
	for (int i = 0; i < nNodes; i++)
	{
		maxBel = -1;
		p_nodeBel = nodeBel + i;
		for (int j = 0; j < nStates[i]; i++)
		{
			if (p_nodeBel[0] > maxBel)
			{
				maxBel = p_nodeBel[0];
				labels[i] = j;
			}
			p_nodeBel += nNodes;
		}

	}

	for (int i = 0; i < nNodes; i++)
		labels[i]++;

	UNPROTECT(6);
	return(_labels);
}
