#include "CRF.h"

/* Tree BP */

void CRF::TreeBP(double *messages_1, double *messages_2, bool maximize)
{
	for (int i = 0; i < maxState * nEdges; i++)
		messages_1[i] = messages_2[i] = 0;
	int *nWaiting = (int *) R_alloc(nNodes, sizeof(int));
	int **waiting = (int **) R_alloc(nNodes, sizeof(int *));
	int *sent = (int *) R_alloc(nNodes, sizeof(int));
	int senderQueueHead, senderQueueTail;
	int *senderQueue = (int *) R_alloc(nNodes * 2, sizeof(int *));
	senderQueueHead = senderQueueTail = 0;
	for (int i = 0; i < nNodes; i++)
	{
		nWaiting[i] = nAdj[i];
		waiting[i] = (int *) R_alloc(nAdj[i], sizeof(int));
		for (int j = 0; j < nAdj[i]; j++)
			waiting[i][j] = 1;
		sent[i] = -1;
		if (nAdj[i] == 1)
			senderQueue[senderQueueTail++] = i;
	}

	int nReceiverQueue;
	int *receiverQueue = (int *) R_alloc(nNodes, sizeof(int *));
	double *incoming = (double *) R_alloc(maxState, sizeof(double));

	int s, r, e, n;
	double mesg, sumMesg, *p_messages;

	while (senderQueueHead < senderQueueTail)
	{
		R_CheckUserInterrupt();

		s = senderQueue[senderQueueHead++];
		if (sent[s] == -2) continue;

		nReceiverQueue = 0;
		if (nWaiting[s] == 1)
		{
			for (int i = 0; i < nAdj[s]; i++)
			{
				if (waiting[s][i])
				{
					receiverQueue[nReceiverQueue++] = i;
					sent[s] = nAdj[s] == 1 ? -2 : i;
					break;
				}
			}
		}
		else
		{
			for (int i = 0; i < nAdj[s]; i++)
				if (sent[s] != i)
					receiverQueue[nReceiverQueue++] = i;
			sent[s] = -2;
		}

		for (int i = 0; i < nReceiverQueue; i++)
		{
			n = receiverQueue[i];
			r = AdjNodes(s, n);

			for (int j = 0; j < nAdj[r]; j++)
				if (AdjNodes(r, j) == s)
				{
					waiting[r][j] = 0;
					nWaiting[r]--;
					break;
				}

			if (sent[r] != -2 && nWaiting[r] <= 1)
				senderQueue[senderQueueTail++] = r;

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
	}
}
