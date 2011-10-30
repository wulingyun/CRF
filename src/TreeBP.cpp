#include "CRF.h"

/* Tree BP */

void CRF::TreeBP(double *messages_1, double *messages_2, bool maximize)
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
	double mesg, sumMesg, *p_messages;

	int done = 0;
	while (!done)
	{
		R_CheckUserInterrupt();

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
					r = AdjNodes(s, n);

					unsent[s][n] = 0;
					nUnsent[s]--;

					for (int j = 0; j < nAdj[r]; j++)
						if (AdjNodes(r, j) == s)
						{
							waiting[r][j] = 0;
							nWaiting[r]--;
							break;
						}

					/* gather incoming messages */

					for (int j = 0; j < nStates[s]; j++)
						incoming[j] = NodePot(s, j);
					for (int j = 0; j < nAdj[s]; j++)
					{
						if (j != n)
						{
							e = AdjEdges(s, j);
							if (EdgesBegin(e) == s)
								p_messages = messages_1;
							else
								p_messages = messages_2;
							p_messages += maxState * e;
							for (int k = 0; k < nStates[s]; k++)
								incoming[k] *= p_messages[k];
						}
					}

					/* send messages */

					e = AdjEdges(s, n);
					sumMesg = 0;
					if (EdgesBegin(e) == s)
					{
						p_messages = messages_2 + maxState * e;
						for (int j = 0; j < nStates[r]; j++)
						{
							p_messages[j] = 0;
							if (maximize)
							{
								for (int k = 0; k < nStates[s]; k++)
								{
									mesg = incoming[k] * EdgePot(e, k, j);
									if (mesg > p_messages[j])
										p_messages[j] = mesg;
								}
							}
							else
							{
								for (int k = 0; k < nStates[s]; k++)
									p_messages[j] += incoming[k] * EdgePot(e, k, j);
							}
							sumMesg += p_messages[j];
						}
					}
					else
					{
						p_messages = messages_1 + maxState * e;
						for (int j = 0; j < nStates[r]; j++)
						{
							p_messages[j] = 0;
							if (maximize)
							{
								for (int k = 0; k < nStates[s]; k++)
								{
									mesg = incoming[k] * EdgePot(e, j, k);
									if (mesg > p_messages[j])
										p_messages[j] = mesg;
								}
							}
							else
							{
								for (int k = 0; k < nStates[s]; k++)
								{
									p_messages[j] += incoming[k] * EdgePot(e, j, k);
								}
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
