#include "CRF.h"

/* Tree BP */

void CRF::TreeBP(bool maximize)
{
	for (int i = 0; i < nEdges; i++)
		for (int j = 0; j < maxState; j++)
		{
			messages[0][i][j] = messages[1][i][j] = 1;
		}

	int *nWaiting = (int *) R_alloc(nNodes, sizeof(int));
	int **waiting = (int **) R_allocArray2<int>(nNodes, nAdj);
	int *sent = (int *) R_alloc(nNodes, sizeof(int));
	int senderHead, senderTail, nReceiver;
	int *sender = (int *) R_alloc(nNodes * 2, sizeof(int));
	int *receiver = (int *) R_alloc(nNodes, sizeof(int));

	senderHead = senderTail = nReceiver = 0;
	for (int i = 0; i < nNodes; i++)
	{
		nWaiting[i] = nAdj[i];
		for (int j = 0; j < nAdj[i]; j++)
			waiting[i][j] = 1;
		sent[i] = -1;
		if (nAdj[i] == 1)
			sender[senderTail++] = i;
		for (int j = 0; j < nStates[i]; j++)
			NodeBel(i, j) = NodePot(i, j);
	}

	int s, r, e, n;
	double m, *msg, sumMsg;
	double *outgoing = (double *) R_alloc(maxState, sizeof(double));

	while (senderHead < senderTail)
	{
		R_CheckUserInterrupt();

		s = sender[senderHead++];
		if (sent[s] == -2) continue;

		nReceiver = 0;
		if (nWaiting[s] == 1)
		{
			for (int i = 0; i < nAdj[s]; i++)
			{
				if (waiting[s][i])
				{
					receiver[nReceiver++] = i;
					sent[s] = nAdj[s] == 1 ? -2 : i;
					break;
				}
			}
		}
		else
		{
			for (int i = 0; i < nAdj[s]; i++)
				if (sent[s] != i)
					receiver[nReceiver++] = i;
			sent[s] = -2;
		}

		for (int i = 0; i < nReceiver; i++)
		{
			n = receiver[i];
			r = AdjNodes(s, n);

			for (int j = 0; j < nAdj[r]; j++)
				if (AdjNodes(r, j) == s)
				{
					waiting[r][j] = 0;
					nWaiting[r]--;
					break;
				}

			if (sent[r] != -2 && nWaiting[r] <= 1)
				sender[senderTail++] = r;

			/* send messages */

			e = AdjEdges(s, n);
			sumMsg = 0;
			if (EdgesBegin(e) == s)
			{
				for (int j = 0; j < nStates[s]; j++)
					outgoing[j] = messages[0][e][j] == 0 ? 0 : NodeBel(s, j) / messages[0][e][j];
				msg = messages[1][e];
				for (int j = 0; j < nStates[r]; j++)
				{
					msg[j] = 0;
					if (maximize)
					{
						for (int k = 0; k < nStates[s]; k++)
						{
							m = outgoing[k] * EdgePot(e, k, j);
							if (m > msg[j])
								msg[j] = m;
						}
					}
					else
					{
						for (int k = 0; k < nStates[s]; k++)
							msg[j] += outgoing[k] * EdgePot(e, k, j);
					}
					sumMsg += msg[j];
				}
			}
			else
			{
				for (int j = 0; j < nStates[s]; j++)
					outgoing[j] = messages[1][e][j] == 0 ? 0 : NodeBel(s, j) / messages[1][e][j];
				msg = messages[0][e];
				for (int j = 0; j < nStates[r]; j++)
				{
					msg[j] = 0;
					if (maximize)
					{
						for (int k = 0; k < nStates[s]; k++)
						{
							m = outgoing[k] * EdgePot(e, j, k);
							if (m > msg[j])
								msg[j] = m;
						}
					}
					else
					{
						for (int k = 0; k < nStates[s]; k++)
						{
							msg[j] += outgoing[k] * EdgePot(e, j, k);
						}
					}
					sumMsg += msg[j];
				}
			}
			for (int j = 0; j < nStates[r]; j++)
			{
				msg[j] /= sumMsg;
				NodeBel(r, j) *= msg[j];
			}
		}
	}

	double sumBel;
	for (int i = 0; i < nNodes; i++)
	{
		sumBel = 0;
		for (int j = 0; j < nStates[i]; j++)
			sumBel += NodeBel(i, j);
		for (int j = 0; j < nStates[i]; j++)
			NodeBel(i, j) /= sumBel;
	}
}
