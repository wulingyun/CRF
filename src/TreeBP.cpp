#include "CRF.h"

/* Tree BP */

void CRF::TreeBP(bool maximize)
{
	messages = (double ***) R_allocArray<double>(2, nEdges, maxState);
	for (int i = 0; i < nEdges; i++)
		for (int j = 0; j < maxState; j++)
			messages[0][i][j] = messages[1][i][j] = 1;

	int *nWaiting = (int *) R_alloc(nNodes, sizeof(int));
	int **waiting = (int **) R_allocArray2<int>(nNodes, nAdj);
	int *sent = (int *) R_alloc(nNodes, sizeof(int));
	int senderHead, senderTail, nReceiver;
	int *sender = (int *) R_alloc(nNodes * 2, sizeof(int));
	int *receiver = (int *) R_alloc(nNodes, sizeof(int));

	double sumBel;
	senderHead = senderTail = nReceiver = 0;
	for (int i = 0; i < nNodes; i++)
	{
		nWaiting[i] = nAdj[i];
		for (int j = 0; j < nAdj[i]; j++)
			waiting[i][j] = 1;
		sent[i] = -1;
		if (nAdj[i] == 1)
			sender[senderTail++] = i;
		sumBel = 0;
		for (int j = 0; j < nStates[i]; j++)
		  sumBel += NodeBel(i, j) = NodePot(i, j);
		for (int j = 0; j < nStates[i]; j++)
		  NodeBel(i, j) /= sumBel;
	}

	int s, r, e, n;
	double *msg;
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

		/* send messages */

		for (int i = 0; i < nReceiver; i++)
		{
			n = receiver[i];
			r = AdjNodes(s, n);
			e = AdjEdges(s, n);

			for (int j = 0; j < nAdj[r]; j++)
				if (AdjNodes(r, j) == s)
				{
					waiting[r][j] = 0;
					nWaiting[r]--;
					break;
				}

			if (sent[r] != -2 && nWaiting[r] <= 1)
				sender[senderTail++] = r;

			if (maximize)
				msg = SendMessagesMax(s, r, e, outgoing, messages);
			else
				msg = SendMessagesSum(s, r, e, outgoing, messages);

			sumBel = 0;
			for (int j = 0; j < nStates[r]; j++)
			  sumBel += NodeBel(r, j) *= msg[j];
			for (int j = 0; j < nStates[r]; j++)
			  NodeBel(r, j) /= sumBel;
		}
	}
}
