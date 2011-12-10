#include "CRF.h"

/* Loopy BP */

void CRF::LoopyBP(int maxIter, double cutoff, int verbose, bool maximize)
{
	int dim[] = {2, nEdges, maxState};
	double ***old_messages = (double ***) allocArray<double, 3>(dim);

	for (int i = 0; i < nEdges; i++)
		for (int j = 0; j < maxState; j++)
		{
			messages[0][i][j] = old_messages[0][i][j] = 0;
			messages[1][i][j] = old_messages[1][i][j] = 0;
		}

	double *incoming = (double *) R_alloc(maxState, sizeof(double));
	double *outgoing = (double *) R_alloc(maxState, sizeof(double));

	int s, r, e, n;
	double m, *msg, sumMsg;

	for (int i = 0; i < nEdges; i++)
	{
		n = EdgesBegin(i);
		for (int j = 0; j < nStates[n]; j++)
			messages[0][i][j] = 1.0 / nStates[n];
		n = EdgesEnd(i);
		for (int j = 0; j < nStates[n]; j++)
			messages[1][i][j] = 1.0 / nStates[n];
	}

	double difference = 0;
	for (int iter = 1; iter <= maxIter; iter++)
	{
		R_CheckUserInterrupt();

		swap(old_messages, messages);

		for (s = 0; s < nNodes; s++)
		{
			/* gather incoming messages */

			for (int i = 0; i < nStates[s]; i++)
				incoming[i] = NodePot(s, i);
			for (int i = 0; i < nAdj[s]; i++)
			{
				e = AdjEdges(s, i);
				if (EdgesBegin(e) == s)
					msg = old_messages[0][e];
				else
					msg = old_messages[1][e];
				for (int k = 0; k < nStates[s]; k++)
					incoming[k] *= msg[k];
			}

			/* send messages */

			for (int i = 0; i < nAdj[s]; i++)
			{
				r = AdjNodes(s, i);
				e = AdjEdges(s, i);

				if (EdgesBegin(e) == s)
					msg = old_messages[0][e];
				else
					msg = old_messages[1][e];
				for (int k = 0; k < nStates[s]; k++)
					outgoing[k] = msg[k] == 0 ? 0 : incoming[k] / msg[k];

				sumMsg = 0;
				if (EdgesBegin(e) == s)
				{
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
								msg[j] += outgoing[k] * EdgePot(e, j, k);
						}
						sumMsg += msg[j];
					}
				}
				for (int j = 0; j < nStates[r]; j++)
					msg[j] /= sumMsg;
			}
		}

		difference = 0;
		for (int i = 0; i < nEdges; i++)
			for (int j = 0; j < maxState; j++)
			{
				difference += fabs(messages[0][i][j] - old_messages[0][i][j]);
				difference += fabs(messages[1][i][j] - old_messages[1][i][j]);
			}
		if (verbose)
			Rprintf("LBP: Iteration %d, Difference = %f\n", iter, difference);
		if (difference <= cutoff)
			break;
	}

	if (difference > cutoff)
		warning("Loopy BP did not converge in %d iterations! (diff = %f)", maxIter, difference);
}
