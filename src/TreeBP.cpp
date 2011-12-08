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
	int nSenderQueue, senderQueueHead, senderQueueTail;
	int *senderQueue = (int *) R_alloc(nNodes * 2, sizeof(int *));
	nSenderQueue = senderQueueHead = senderQueueTail = 0;
	for (int i = 0; i < nNodes; i++)
	{
		nWaiting[i] = nUnsent[i] = nAdj[i];
		waiting[i] = (int *) R_alloc(nAdj[i], sizeof(int));
		unsent[i] = (int *) R_alloc(nAdj[i], sizeof(int));
		for (int j = 0; j < nAdj[i]; j++)
			waiting[i][j] = unsent[i][j] = 1;
		if (nAdj[i] == 1)
		{
			nSenderQueue++;
			senderQueue[senderQueueTail++] = i;
		}
	}

	int nReceiverQueue;
	int *receiverQueue = (int *) R_alloc(nNodes, sizeof(int *));
	double *incoming = (double *) R_alloc(maxState, sizeof(double));

	int s, r, e, n;
	double mesg, sumMesg, *p_messages;

	while (nSenderQueue > 0)
	{
		R_CheckUserInterrupt();

		nSenderQueue--;
		s = senderQueue[senderQueueHead++];

		nReceiverQueue = 0;
		if (nWaiting[s] == 1)
		{
			for (int i = 0; i < nAdj[s]; i++)
				if (waiting[s][i] && unsent[s][i])
					receiverQueue[nReceiverQueue++] = i;
		}
		else
		{
			for (int i = 0; i < nAdj[s]; i++)
				if (unsent[s][i])
					receiverQueue[nReceiverQueue++] = i;
		}

		for (int i = 0; i < nReceiverQueue; i++)
		{
			n = receiverQueue[i];
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

			if (nUnsent[r] > 0 && nWaiting[r] <= 1)
			{
				nSenderQueue++;
				senderQueue[senderQueueTail++] = r;
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
	}
}
