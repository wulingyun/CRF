#include "CRF.h"

SEXP Decode_LBP(SEXP _crf, SEXP _maxIter, SEXP _cutoff, SEXP _verbose)
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

	PROTECT(_maxIter = AS_INTEGER(_maxIter));
	int maxIter = INTEGER_POINTER(_maxIter)[0];
	PROTECT(_cutoff = AS_NUMERIC(_cutoff));
	double cutoff = NUMERIC_POINTER(_cutoff)[0];
	PROTECT(_verbose = AS_INTEGER(_verbose));
	int verbose = INTEGER_POINTER(_verbose)[0];

	SEXP _labels;
	PROTECT(_labels = NEW_INTEGER(nNodes));
	int *labels = INTEGER_POINTER(_labels);
	setValues(_labels, labels, -1);

	/* Loopy BP */

	double *messages_1 = (double *) R_alloc(maxState * nEdges, sizeof(double)); // Messages from n2 to n1 at edge (n1, n2)
	double *messages_2 = (double *) R_alloc(maxState * nEdges, sizeof(double)); // Messages from n1 to n2 at edge (n1, n2)
	double *old_messages_1 = (double *) R_alloc(maxState * nEdges, sizeof(double));
	double *old_messages_2 = (double *) R_alloc(maxState * nEdges, sizeof(double));
	for (int i = 0; i < maxState * nEdges; i++)
		messages_1[i] = messages_2[i] = old_messages_1[i] = old_messages_2[i] = 0;

	double *incoming = (double *) R_alloc(maxState, sizeof(double));

	int s, r, e, n;
	double mesg, sumMesg, *p_nodePot, *p_edgePot, *p0_edgePot, *p_messages;

	for (int i = 0; i < nEdges; i++)
	{
		p_messages = messages_1 + maxState * i;
		n = edges[i] - 1;
		for (int j = 0; j < nStates[n]; j++)
			p_messages[j] = 1.0 / nStates[n];
		p_messages = messages_2 + maxState * i;
		n = edges[i + nEdges] - 1;
		for (int j = 0; j < nStates[n]; j++)
			p_messages[j] = 1.0 / nStates[n];
	}

	double difference = 0;
	for (int iter = 1; iter <= maxIter; iter++)
	{
		p_messages = old_messages_1;
		old_messages_1 = messages_1;
		messages_1 = p_messages;
		p_messages = old_messages_2;
		old_messages_2 = messages_2;
		messages_2 = p_messages;

		for (s = 0; s < nNodes; s++)
		{
			for (int i = 0; i < nAdj[s]; i++)
			{
				r = adjNodes[s][i] - 1;

				/* gather incoming messages */

				p_nodePot = nodePot + s;
				for (int j = 0; j < nStates[s]; j++)
				{
					incoming[j] = p_nodePot[0];
					p_nodePot += nNodes;
				}
				for (int j = 0; j < nAdj[s]; j++)
				{
					if (j != i)
					{
						e = adjEdges[s][j] - 1;
						if (edges[e] - 1 == s)
							p_messages = old_messages_1;
						else
							p_messages = old_messages_2;
						p_messages += maxState * e;
						for (int k = 0; k < nStates[s]; k++)
							incoming[k] *= p_messages[k];
					}
				}

				/* send messages */

				e = adjEdges[s][i] - 1;
				sumMesg = 0;
				p0_edgePot = edgePot + maxState * maxState * e;
				if (edges[e] - 1 == s)
				{
					p_messages = messages_2 + maxState * e;
					for (int j = 0; j < nStates[r]; j++)
					{
						p_edgePot = p0_edgePot;
						p0_edgePot += maxState;
						p_messages[j] = -1;
						for (int k = 0; k < nStates[s]; k++)
						{
							mesg = incoming[k] * p_edgePot[k];
							if (mesg > p_messages[j])
								p_messages[j] = mesg;
						}
						sumMesg += p_messages[j];
					}
				}
				else
				{
					p_messages = messages_1 + maxState * e;
					for (int j = 0; j < nStates[r]; j++)
					{
						p_edgePot = p0_edgePot++;
						p_messages[j] = -1;
						for (int k = 0; k < nStates[s]; k++)
						{
							mesg = incoming[k] * p_edgePot[0];
							if (mesg > p_messages[j])
								p_messages[j] = mesg;
							p_edgePot += maxState;
						}
						sumMesg += p_messages[j];
					}
				}
				for (int j = 0; j < nStates[r]; j++)
					p_messages[j] /= sumMesg;
			}
		}

		difference = 0;
		for (int i = 0; i < maxState * nEdges; i++)
		{
			difference += fabs(messages_1[i] - old_messages_1[i]);
			difference += fabs(messages_2[i] - old_messages_2[i]);
		}
		if (verbose)
			Rprintf("LBP: Iteration %d, Difference = %f\n", iter, difference);
		if (difference <= cutoff)
			break;
	}

	if (difference > cutoff)
		warning("Loopy BP did not converge in %d iterations! (diff = %f)", maxIter, difference);

	/* Node beliefs */

	double *nodeBel = (double *) R_alloc(nNodes * maxState, sizeof(double));
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

	/* Labels */

	double maxBel;
	for (int i = 0; i < nNodes; i++)
	{
		maxBel = -1;
		p_nodeBel = nodeBel + i;
		for (int j = 0; j < nStates[i]; j++)
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

	UNPROTECT(14 + nNodes * 2);

	return(_labels);
}
