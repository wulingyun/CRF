#include "CRF.h"

SEXP Infer_Tree(SEXP _crf)
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

	SEXP _nAdj, _adjNodes, _adjEdges, _temp;
	PROTECT(_nAdj = AS_INTEGER(getListElement(_crf, "n.adj")));
	PROTECT(_adjNodes = AS_LIST(getListElement(_crf, "adj.nodes")));
	PROTECT(_adjEdges = AS_LIST(getListElement(_crf, "adj.edges")));
	int *nAdj = INTEGER_POINTER(_nAdj);
	int **adjNodes = (int **) R_alloc(nNodes, sizeof(int *));
	int **adjEdges = (int **) R_alloc(nNodes, sizeof(int *));
	for (int i = 0; i < nNodes; i++)
	{
		PROTECT(_temp = AS_INTEGER(VECTOR_ELT(_adjNodes, i)));
		adjNodes[i] = INTEGER_POINTER(_temp);
		PROTECT(_temp = AS_INTEGER(VECTOR_ELT(_adjEdges, i)));
		adjEdges[i] = INTEGER_POINTER(_temp);
	}

	SEXP _nodePot, _edgePot;
	PROTECT(_nodePot = AS_NUMERIC(getListElement(_crf, "node.pot")));
	PROTECT(_edgePot = AS_NUMERIC(getListElement(_crf, "edge.pot")));
	double *nodePot = NUMERIC_POINTER(_nodePot);
	double *edgePot = NUMERIC_POINTER(_edgePot);

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

	/* Tree BP */

	double *messages_1 = (double *) R_alloc(maxState * nEdges, sizeof(double)); // Messages from n2 to n1 at edge (n1, n2)
	double *messages_2 = (double *) R_alloc(maxState * nEdges, sizeof(double)); // Messages from n1 to n2 at edge (n1, n2)
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
	double *incoming = (double *) R_alloc(maxState, sizeof(double));

	int s, r, e, n;
	double sumMesg, *p_nodePot, *p_edgePot, *p0_edgePot, *p_messages;

	int done = 0;
	while (!done)
	{
		done = 1;
		for (s = 0; s < nNodes; s++)
		{
			if (nUnsent[s] == 0 || nWaiting[s] > 1)
				continue;

			nQueue = 0;
			if (nWaiting[s] == 1)
			{
				for (int i = 0; i < nAdj[s]; i++)
					if (waiting[s][i] && unsent[s][i])
						queue[nQueue++] = i;
			}
			else
			{
				for (int i = 0; i < nAdj[s]; i++)
					if (unsent[s][i])
						queue[nQueue++] = i;
			}

			if (nQueue > 0)
			{
				for (int i = 0; i < nQueue; i++)
				{
					n = queue[i];
					r = adjNodes[s][n] - 1;

					unsent[s][n] = 0;
					nUnsent[s]--;

					for (int j = 0; j < nAdj[r]; j++)
						if (adjNodes[r][j] - 1 == s)
						{
							waiting[r][j] = 0;
							nWaiting[r]--;
							break;
						}

					/* gather incoming messages */

					p_nodePot = nodePot + s;
					for (int j = 0; j < nStates[s]; j++)
					{
						incoming[j] = p_nodePot[0];
						p_nodePot += nNodes;
					}
					for (int j = 0; j < nAdj[s]; j++)
					{
						if (j != n)
						{
							e = adjEdges[s][j] - 1;
							if (edges[e] - 1 == s)
								p_messages = messages_1;
							else
								p_messages = messages_2;
							p_messages += maxState * e;
							for (int k = 0; k < nStates[s]; k++)
								incoming[k] *= p_messages[k];
						}
					}

					/* send messages */

					e = adjEdges[s][n] - 1;
					sumMesg = 0;
					p0_edgePot = edgePot + maxState * maxState * e;
					if (edges[e] - 1 == s)
					{
						p_messages = messages_2 + maxState * e;
						for (int j = 0; j < nStates[r]; j++)
						{
							p_edgePot = p0_edgePot;
							p0_edgePot += maxState;
							p_messages[j] = 0;
							for (int k = 0; k < nStates[s]; k++)
								p_messages[j] += incoming[k] * p_edgePot[k];
							sumMesg += p_messages[j];
						}
					}
					else
					{
						p_messages = messages_1 + maxState * e;
						for (int j = 0; j < nStates[r]; j++)
						{
							p_edgePot = p0_edgePot++;
							p_messages[j] = 0;
							for (int k = 0; k < nStates[s]; k++)
							{
								p_messages[j] += incoming[k] * p_edgePot[0];
								p_edgePot += maxState;
							}
							sumMesg += p_messages[j];
						}
					}
					for (int j = 0; j < nStates[r]; j++)
						p_messages[j] /= sumMesg;
				}
				done = 0;
			}
		}
	}

	/* Node beliefs */

	for (int i = 0; i < length(_nodePot); i++)
		nodeBel[i] = nodePot[i];

	int n1, n2;
	double *p_nodeBel;
	double *p1_messages = messages_1;
	double *p2_messages = messages_2;
	for (int i = 0; i < nEdges; i++)
	{
		n1 = edges[i] - 1;
		n2 = edges[i + nEdges] - 1;
		p_nodeBel = nodeBel + n1;
		for (int j = 0; j < nStates[n1]; j++)
		{
			p_nodeBel[0] *= p1_messages[j];
			p_nodeBel += nNodes;
		}
		p_nodeBel = nodeBel + n2;
		for (int j = 0; j < nStates[n2]; j++)
		{
			p_nodeBel[0] *= p2_messages[j];
			p_nodeBel += nNodes;
		}
		p1_messages += maxState;
		p2_messages += maxState;
	}

	double sumBel;
	for (int i = 0; i < nNodes; i++)
	{
		sumBel = 0;
		p_nodeBel = nodeBel + i;
		for (int j = 0; j < nStates[i]; j++)
		{
			sumBel += p_nodeBel[0];
			p_nodeBel += nNodes;
		}
		p_nodeBel = nodeBel + i;
		for (int j = 0; j < nStates[i]; j++)
		{
			p_nodeBel[0] /= sumBel;
			p_nodeBel += nNodes;
		}
	}

	/* Edge beliefs */

	for (int i = 0; i < length(_edgePot); i++)
		edgeBel[i] = edgePot[i];

	double bel, *p_edgeBel, *p0_edgeBel;
	p0_edgeBel = edgeBel;
	p1_messages = messages_1;
	p2_messages = messages_2;
	for (int i = 0; i < nEdges; i++)
	{
		n1 = edges[i] - 1;
		n2 = edges[i + nEdges] - 1;
		p_nodeBel = nodeBel + n1;
		for (int j = 0; j < nStates[n1]; j++)
		{
			if (p1_messages[j] != 0)
				bel = p_nodeBel[0] / p1_messages[j];
			else
				bel = 0;
			p_edgeBel = p0_edgeBel + j;
			for (int k = 0; k < nStates[n2]; k++)
			{
				p_edgeBel[0] *= bel;
				p_edgeBel += maxState;
			}
			p_nodeBel += nNodes;
		}
		p_nodeBel = nodeBel + n2;
		p_edgeBel = p0_edgeBel;
		for (int j = 0; j < nStates[n2]; j++)
		{
			if (p2_messages[j] != 0)
				bel = p_nodeBel[0] / p2_messages[j];
			else
				bel = 0;
			for (int k = 0; k < nStates[n2]; k++)
				p_edgeBel[k] *= bel;
			p_nodeBel += nNodes;
			p_edgeBel += maxState;
		}

		sumBel = 0;
		p_edgeBel = p0_edgeBel;
		for (int j = 0; j < nStates[n2]; j++)
		{
			for (int k = 0; k < nStates[n1]; k++)
				sumBel += p_edgeBel[k];
			p_edgeBel += maxState;
		}
		p_edgeBel = p0_edgeBel;
		for (int j = 0; j < nStates[n2]; j++)
		{
			for (int k = 0; k < nStates[n1]; k++)
				p_edgeBel[k] /= sumBel;
			p_edgeBel += maxState;
		}

		p0_edgeBel += maxState * maxState;
		p1_messages += maxState;
		p2_messages += maxState;
	}

	/* Bethe free energy */

	double nodeEnergy, nodeEntropy, edgeEnergy, edgeEntropy;
	nodeEnergy = nodeEntropy = edgeEnergy = edgeEntropy = 0;

	double entropy;
	for (int i = 0; i < nNodes; i++)
	{
		entropy = 0;
		p_nodeBel = nodeBel + i;
		p_nodePot = nodePot + i;
		for (int j = 0; j < nStates[i]; j++)
		{
			if (p_nodeBel[0] > 0)
			{
				nodeEnergy -= p_nodeBel[0] * log(p_nodePot[0]);
				entropy += p_nodeBel[0] * log(p_nodeBel[0]);
			}
			p_nodeBel += nNodes;
			p_nodePot += nNodes;
		}
		nodeEntropy += (nAdj[i] - 1) * entropy;
	}

	p0_edgeBel = edgeBel;
	p0_edgePot = edgePot;
	for (int i = 0; i < nEdges; i++)
	{
		n1 = edges[i] - 1;
		n2 = edges[i + nEdges] - 1;
		p_edgeBel = p0_edgeBel;
		p_edgePot = p0_edgePot;
		for (int j = 0; j < nStates[n2]; j++)
		{
			for (int k = 0; k < nStates[n1]; k++)
			{
				if (p_edgeBel[k] > 0)
				{
					edgeEnergy -= p_edgeBel[k] * log(p_edgePot[k]);
					edgeEntropy -= p_edgeBel[k] * log(p_edgeBel[k]);
				}
			}
			p_edgeBel += maxState;
			p_edgePot += maxState;
		}
		p0_edgeBel += maxState * maxState;
		p0_edgePot += maxState * maxState;
	}

	*logZ = - nodeEnergy + nodeEntropy - edgeEnergy + edgeEntropy;

	SEXP _belief;
	PROTECT(_belief = NEW_LIST(3));
	setListElement(_belief, 0, "node.bel", _nodeBel);
	setListElement(_belief, 1, "edge.bel", _edgeBel);
	setListElement(_belief, 2, "logZ", _logZ);

	UNPROTECT(14 + nNodes * 2);

	return(_belief);
}
