#include "CRF.h"

/* Tree BP (sum product) */

void CRF::TreeBP(double *messages_1, double *messages_2)
{
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
}

/* Tree BP (max product) */

void CRF::TreeBP_max(double *messages_1, double *messages_2)
{
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
	double mesg, sumMesg, *p_nodePot, *p_edgePot, *p0_edgePot, *p_messages;

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
				done = 0;
			}
		}
	}
}

/* Loopy BP (sum product) */

void CRF::LoopyBP(double *messages_1, double *messages_2, int maxIter, double cutoff, int verbose)
{
	double *old_messages_1 = (double *) R_alloc(maxState * nEdges, sizeof(double));
	double *old_messages_2 = (double *) R_alloc(maxState * nEdges, sizeof(double));
	for (int i = 0; i < maxState * nEdges; i++)
		messages_1[i] = messages_2[i] = old_messages_1[i] = old_messages_2[i] = 0;

	double *incoming = (double *) R_alloc(maxState, sizeof(double));

	int s, r, e, n;
	double sumMesg, *p_nodePot, *p_edgePot, *p0_edgePot, *p_messages;

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
}

/* Loopy BP (max product) */

void CRF::LoopyBP_max(double *messages_1, double *messages_2, int maxIter, double cutoff, int verbose)
{
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
}

/* Node beliefs */

void CRF::Message2NodeBelief(double *messages_1, double *messages_2)
{
	for (int i = 0; i < length(_nodePot); i++)
		nodeBel[i] = nodePot[i];

	int n1, n2;
	double sumBel, *p_nodeBel;
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
}

/* Edge beliefs */

void CRF::Message2EdgeBelief(double *messages_1, double *messages_2)
{
	for (int i = 0; i < length(_edgePot); i++)
		edgeBel[i] = edgePot[i];

	int n1, n2;
	double bel, sumBel, *p_nodeBel, *p_edgeBel;
	double *p0_edgeBel = edgeBel;
	double *p1_messages = messages_1;
	double *p2_messages = messages_2;
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
}

/* Decoding by max of marginals */

void CRF::MaxOfMarginals()
{
	double maxBel, *p_nodeBel;
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
}

/* Bethe free energy */

void CRF::BetheFreeEnergy()
{
	double nodeEnergy, nodeEntropy, edgeEnergy, edgeEntropy;
	nodeEnergy = nodeEntropy = edgeEnergy = edgeEntropy = 0;

	double entropy;
	double *p_nodeBel, *p_nodePot;
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

	int n1, n2;
	double *p_edgeBel, *p0_edgeBel = edgeBel;
	double *p_edgePot, *p0_edgePot = edgePot;
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
}
