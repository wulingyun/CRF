#include "CRF.h"

/* Loopy BP */

void CRF::LoopyBP(double *messages_1, double *messages_2, int maxIter, double cutoff, int verbose, bool maximize)
{
	double *old_messages_1 = (double *) R_alloc(maxState * nEdges, sizeof(double));
	double *old_messages_2 = (double *) R_alloc(maxState * nEdges, sizeof(double));
	for (int i = 0; i < maxState * nEdges; i++)
		messages_1[i] = messages_2[i] = old_messages_1[i] = old_messages_2[i] = 0;

	double *incoming = (double *) R_alloc(maxState, sizeof(double));
	double *outgoing = (double *) R_alloc(maxState, sizeof(double));

	int s, r, e, n;
	double mesg, sumMesg, *p_messages;

	for (int i = 0; i < nEdges; i++)
	{
		p_messages = messages_1 + maxState * i;
		n = EdgesBegin(i);
		for (int j = 0; j < nStates[n]; j++)
			p_messages[j] = 1.0 / nStates[n];
		p_messages = messages_2 + maxState * i;
		n = EdgesEnd(i);
		for (int j = 0; j < nStates[n]; j++)
			p_messages[j] = 1.0 / nStates[n];
	}

	double difference = 0;
	for (int iter = 1; iter <= maxIter; iter++)
	{
		R_CheckUserInterrupt();

		p_messages = old_messages_1;
		old_messages_1 = messages_1;
		messages_1 = p_messages;
		p_messages = old_messages_2;
		old_messages_2 = messages_2;
		messages_2 = p_messages;

		for (s = 0; s < nNodes; s++)
		{
			/* gather incoming messages */

			for (int i = 0; i < nStates[s]; i++)
				incoming[i] = NodePot(s, i);
			for (int i = 0; i < nAdj[s]; i++)
			{
				e = AdjEdges(s, i);
				if (EdgesBegin(e) == s)
					p_messages = old_messages_1;
				else
					p_messages = old_messages_2;
				p_messages += maxState * e;
				for (int k = 0; k < nStates[s]; k++)
					incoming[k] *= p_messages[k];
			}

			/* send messages */

			for (int i = 0; i < nAdj[s]; i++)
			{
				r = AdjNodes(s, i);
				e = AdjEdges(s, i);

				if (EdgesBegin(e) == s)
					p_messages = old_messages_1;
				else
					p_messages = old_messages_2;
				p_messages += maxState * e;
				for (int k = 0; k < nStates[s]; k++)
					outgoing[k] = p_messages[k] == 0 ? 0 : incoming[k] / p_messages[k];

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
								mesg = outgoing[k] * EdgePot(e, k, j);
								if (mesg > p_messages[j])
									p_messages[j] = mesg;
							}
						}
						else
						{
							for (int k = 0; k < nStates[s]; k++)
								p_messages[j] += outgoing[k] * EdgePot(e, k, j);
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
								mesg = outgoing[k] * EdgePot(e, j, k);
								if (mesg > p_messages[j])
									p_messages[j] = mesg;
							}
						}
						else
						{
							for (int k = 0; k < nStates[s]; k++)
								p_messages[j] += outgoing[k] * EdgePot(e, j, k);
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
