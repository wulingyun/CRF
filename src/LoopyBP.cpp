#include "CRF.h"

/* Loopy BP */

void CRF::LoopyBP(int maxIter, double cutoff, int verbose, bool maximize)
{
	messages = (double ***) R_allocArray<double>(2, nEdges, maxState);
	double ***old_messages = (double ***) R_allocArray<double>(2, nEdges, maxState);

	for (int i = 0; i < nEdges; i++)
		for (int j = 0; j < maxState; j++)
		{
			messages[0][i][j] = old_messages[0][i][j] = 0;
			messages[1][i][j] = old_messages[1][i][j] = 0;
		}

	double *outgoing = (double *) R_alloc(maxState, sizeof(double));

	int s, r, e, n;

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

		swap2(old_messages, messages);

		for (s = 0; s < nNodes; s++)
		{
      GatherIncomingMessages(s, old_messages);

			/* send messages */

			for (int i = 0; i < nAdj[s]; i++)
			{
				r = AdjNodes(s, i);
				e = AdjEdges(s, i);

				if (maximize)
					ComputeMessagesMax(s, r, e, outgoing, old_messages, messages);
				else
					ComputeMessagesSum(s, r, e, outgoing, old_messages, messages);
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
