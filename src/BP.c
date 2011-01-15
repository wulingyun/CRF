#include "CRF.h"

void TreeBP(CRFinfo *crf, double *messages_1, double *messages_2)
{
	for (int i = 0; i < crf->maxState * crf->nEdges; i++)
		messages_1[i] = messages_2[i] = 0;
	int *nWaiting = (int *) R_alloc(crf->nNodes, sizeof(int));
	int **waiting = (int **) R_alloc(crf->nNodes, sizeof(int *));
	int *nUnsent = (int *) R_alloc(crf->nNodes, sizeof(int));
	int **unsent = (int **) R_alloc(crf->nNodes, sizeof(int *));
	for (int i = 0; i < crf->nNodes; i++)
	{
		nWaiting[i] = nUnsent[i] = crf->nAdj[i];
		waiting[i] = (int *) R_alloc(crf->nAdj[i], sizeof(int));
		unsent[i] = (int *) R_alloc(crf->nAdj[i], sizeof(int));
		for (int j = 0; j < crf->nAdj[i]; j++)
			waiting[i][j] = unsent[i][j] = 1;
	}

	int nQueue;
	int *queue = (int *) R_alloc(crf->nNodes, sizeof(int *));
	double *incoming = (double *) R_alloc(crf->maxState, sizeof(double));

	int s, r, e, n;
	double sumMesg, *p_nodePot, *p_edgePot, *p0_edgePot, *p_messages;

	int done = 0;
	while (!done)
	{
		done = 1;
		for (s = 0; s < crf->nNodes; s++)
		{
			if (nUnsent[s] == 0 || nWaiting[s] > 1)
				continue;

			nQueue = 0;
			if (nWaiting[s] == 1)
			{
				for (int i = 0; i < crf->nAdj[s]; i++)
					if (waiting[s][i] && unsent[s][i])
						queue[nQueue++] = i;
			}
			else
			{
				for (int i = 0; i < crf->nAdj[s]; i++)
					if (unsent[s][i])
						queue[nQueue++] = i;
			}

			if (nQueue > 0)
			{
				for (int i = 0; i < nQueue; i++)
				{
					n = queue[i];
					r = crf->adjNodes[s][n] - 1;

					unsent[s][n] = 0;
					nUnsent[s]--;

					for (int j = 0; j < crf->nAdj[r]; j++)
						if (crf->adjNodes[r][j] - 1 == s)
						{
							waiting[r][j] = 0;
							nWaiting[r]--;
							break;
						}

					/* gather incoming messages */

					p_nodePot = crf->nodePot + s;
					for (int j = 0; j < crf->nStates[s]; j++)
					{
						incoming[j] = p_nodePot[0];
						p_nodePot += crf->nNodes;
					}
					for (int j = 0; j < crf->nAdj[s]; j++)
					{
						if (j != n)
						{
							e = crf->adjEdges[s][j] - 1;
							if (crf->edges[e] - 1 == s)
								p_messages = messages_1;
							else
								p_messages = messages_2;
							p_messages += crf->maxState * e;
							for (int k = 0; k < crf->nStates[s]; k++)
								incoming[k] *= p_messages[k];
						}
					}

					/* send messages */

					e = crf->adjEdges[s][n] - 1;
					sumMesg = 0;
					p0_edgePot = crf->edgePot + crf->maxState * crf->maxState * e;
					if (crf->edges[e] - 1 == s)
					{
						p_messages = messages_2 + crf->maxState * e;
						for (int j = 0; j < crf->nStates[r]; j++)
						{
							p_edgePot = p0_edgePot;
							p0_edgePot += crf->maxState;
							p_messages[j] = 0;
							for (int k = 0; k < crf->nStates[s]; k++)
								p_messages[j] += incoming[k] * p_edgePot[k];
							sumMesg += p_messages[j];
						}
					}
					else
					{
						p_messages = messages_1 + crf->maxState * e;
						for (int j = 0; j < crf->nStates[r]; j++)
						{
							p_edgePot = p0_edgePot++;
							p_messages[j] = 0;
							for (int k = 0; k < crf->nStates[s]; k++)
							{
								p_messages[j] += incoming[k] * p_edgePot[0];
								p_edgePot += crf->maxState;
							}
							sumMesg += p_messages[j];
						}
					}
					for (int j = 0; j < crf->nStates[r]; j++)
						p_messages[j] /= sumMesg;
				}
				done = 0;
			}
		}
	}
}

void TreeBP_max(CRFinfo *crf, double *messages_1, double *messages_2)
{
	for (int i = 0; i < crf->maxState * crf->nEdges; i++)
		messages_1[i] = messages_2[i] = 0;
	int *nWaiting = (int *) R_alloc(crf->nNodes, sizeof(int));
	int **waiting = (int **) R_alloc(crf->nNodes, sizeof(int *));
	int *nUnsent = (int *) R_alloc(crf->nNodes, sizeof(int));
	int **unsent = (int **) R_alloc(crf->nNodes, sizeof(int *));
	for (int i = 0; i < crf->nNodes; i++)
	{
		nWaiting[i] = nUnsent[i] = crf->nAdj[i];
		waiting[i] = (int *) R_alloc(crf->nAdj[i], sizeof(int));
		unsent[i] = (int *) R_alloc(crf->nAdj[i], sizeof(int));
		for (int j = 0; j < crf->nAdj[i]; j++)
			waiting[i][j] = unsent[i][j] = 1;
	}

	int nQueue;
	int *queue = (int *) R_alloc(crf->nNodes, sizeof(int *));
	double *incoming = (double *) R_alloc(crf->maxState, sizeof(double));

	int s, r, e, n;
	double mesg, sumMesg, *p_nodePot, *p_edgePot, *p0_edgePot, *p_messages;

	int done = 0;
	while (!done)
	{
		done = 1;
		for (s = 0; s < crf->nNodes; s++)
		{
			if (nUnsent[s] == 0 || nWaiting[s] > 1)
				continue;

			nQueue = 0;
			if (nWaiting[s] == 1)
			{
				for (int i = 0; i < crf->nAdj[s]; i++)
					if (waiting[s][i] && unsent[s][i])
						queue[nQueue++] = i;
			}
			else
			{
				for (int i = 0; i < crf->nAdj[s]; i++)
					if (unsent[s][i])
						queue[nQueue++] = i;
			}

			if (nQueue > 0)
			{
				for (int i = 0; i < nQueue; i++)
				{
					n = queue[i];
					r = crf->adjNodes[s][n] - 1;

					unsent[s][n] = 0;
					nUnsent[s]--;

					for (int j = 0; j < crf->nAdj[r]; j++)
						if (crf->adjNodes[r][j] - 1 == s)
						{
							waiting[r][j] = 0;
							nWaiting[r]--;
							break;
						}

					/* gather incoming messages */

					p_nodePot = crf->nodePot + s;
					for (int j = 0; j < crf->nStates[s]; j++)
					{
						incoming[j] = p_nodePot[0];
						p_nodePot += crf->nNodes;
					}
					for (int j = 0; j < crf->nAdj[s]; j++)
					{
						if (j != n)
						{
							e = crf->adjEdges[s][j] - 1;
							if (crf->edges[e] - 1 == s)
								p_messages = messages_1;
							else
								p_messages = messages_2;
							p_messages += crf->maxState * e;
							for (int k = 0; k < crf->nStates[s]; k++)
								incoming[k] *= p_messages[k];
						}
					}

					/* send messages */

					e = crf->adjEdges[s][n] - 1;
					sumMesg = 0;
					p0_edgePot = crf->edgePot + crf->maxState * crf->maxState * e;
					if (crf->edges[e] - 1 == s)
					{
						p_messages = messages_2 + crf->maxState * e;
						for (int j = 0; j < crf->nStates[r]; j++)
						{
							p_edgePot = p0_edgePot;
							p0_edgePot += crf->maxState;
							p_messages[j] = -1;
							for (int k = 0; k < crf->nStates[s]; k++)
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
						p_messages = messages_1 + crf->maxState * e;
						for (int j = 0; j < crf->nStates[r]; j++)
						{
							p_edgePot = p0_edgePot++;
							p_messages[j] = -1;
							for (int k = 0; k < crf->nStates[s]; k++)
							{
								mesg = incoming[k] * p_edgePot[0];
								if (mesg > p_messages[j])
									p_messages[j] = mesg;
								p_edgePot += crf->maxState;
							}
							sumMesg += p_messages[j];
						}
					}
					for (int j = 0; j < crf->nStates[r]; j++)
						p_messages[j] /= sumMesg;
				}
				done = 0;
			}
		}
	}
}

/* Node beliefs */

void Message2NodeBelief(CRFinfo *crf, double *messages_1, double *messages_2, double *nodeBel)
{
	for (int i = 0; i < length(crf->_nodePot); i++)
		nodeBel[i] = crf->nodePot[i];

	int n1, n2;
	double sumBel, *p_nodeBel;
	double *p1_messages = messages_1;
	double *p2_messages = messages_2;
	for (int i = 0; i < crf->nEdges; i++)
	{
		n1 = crf->edges[i] - 1;
		n2 = crf->edges[i + crf->nEdges] - 1;
		p_nodeBel = nodeBel + n1;
		for (int j = 0; j < crf->nStates[n1]; j++)
		{
			p_nodeBel[0] *= p1_messages[j];
			p_nodeBel += crf->nNodes;
		}
		p_nodeBel = nodeBel + n2;
		for (int j = 0; j < crf->nStates[n2]; j++)
		{
			p_nodeBel[0] *= p2_messages[j];
			p_nodeBel += crf->nNodes;
		}
		p1_messages += crf->maxState;
		p2_messages += crf->maxState;
	}

	for (int i = 0; i < crf->nNodes; i++)
	{
		sumBel = 0;
		p_nodeBel = nodeBel + i;
		for (int j = 0; j < crf->nStates[i]; j++)
		{
			sumBel += p_nodeBel[0];
			p_nodeBel += crf->nNodes;
		}
		p_nodeBel = nodeBel + i;
		for (int j = 0; j < crf->nStates[i]; j++)
		{
			p_nodeBel[0] /= sumBel;
			p_nodeBel += crf->nNodes;
		}
	}
}

/* Edge beliefs */

void Message2EdgeBelief(CRFinfo *crf, double *messages_1, double *messages_2, double *nodeBel, double *edgeBel)
{
	for (int i = 0; i < length(crf->_edgePot); i++)
		edgeBel[i] = crf->edgePot[i];

	int n1, n2;
	double bel, sumBel, *p_nodeBel, *p_edgeBel;
	double *p0_edgeBel = edgeBel;
	double *p1_messages = messages_1;
	double *p2_messages = messages_2;
	for (int i = 0; i < crf->nEdges; i++)
	{
		n1 = crf->edges[i] - 1;
		n2 = crf->edges[i + crf->nEdges] - 1;
		p_nodeBel = nodeBel + n1;
		for (int j = 0; j < crf->nStates[n1]; j++)
		{
			if (p1_messages[j] != 0)
				bel = p_nodeBel[0] / p1_messages[j];
			else
				bel = 0;
			p_edgeBel = p0_edgeBel + j;
			for (int k = 0; k < crf->nStates[n2]; k++)
			{
				p_edgeBel[0] *= bel;
				p_edgeBel += crf->maxState;
			}
			p_nodeBel += crf->nNodes;
		}
		p_nodeBel = nodeBel + n2;
		p_edgeBel = p0_edgeBel;
		for (int j = 0; j < crf->nStates[n2]; j++)
		{
			if (p2_messages[j] != 0)
				bel = p_nodeBel[0] / p2_messages[j];
			else
				bel = 0;
			for (int k = 0; k < crf->nStates[n2]; k++)
				p_edgeBel[k] *= bel;
			p_nodeBel += crf->nNodes;
			p_edgeBel += crf->maxState;
		}

		sumBel = 0;
		p_edgeBel = p0_edgeBel;
		for (int j = 0; j < crf->nStates[n2]; j++)
		{
			for (int k = 0; k < crf->nStates[n1]; k++)
				sumBel += p_edgeBel[k];
			p_edgeBel += crf->maxState;
		}
		p_edgeBel = p0_edgeBel;
		for (int j = 0; j < crf->nStates[n2]; j++)
		{
			for (int k = 0; k < crf->nStates[n1]; k++)
				p_edgeBel[k] /= sumBel;
			p_edgeBel += crf->maxState;
		}

		p0_edgeBel += crf->maxState * crf->maxState;
		p1_messages += crf->maxState;
		p2_messages += crf->maxState;
	}
}