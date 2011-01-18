#ifdef CRF_MAX

/* Tree BP */

#if CRF_MAX == 1
void CRF::TreeBP_max(double *messages_1, double *messages_2)
#else
void CRF::TreeBP(double *messages_1, double *messages_2)
#endif
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
#if CRF_MAX == 1
	double mesg;
#endif

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
#if CRF_MAX == 1
							p_messages[j] = -1;
							for (int k = 0; k < nStates[s]; k++)
							{
								mesg = incoming[k] * p_edgePot[k];
								if (mesg > p_messages[j])
									p_messages[j] = mesg;
							}
#else
							p_messages[j] = 0;
							for (int k = 0; k < nStates[s]; k++)
								p_messages[j] += incoming[k] * p_edgePot[k];
#endif
							sumMesg += p_messages[j];
						}
					}
					else
					{
						p_messages = messages_1 + maxState * e;
						for (int j = 0; j < nStates[r]; j++)
						{
							p_edgePot = p0_edgePot++;
#if CRF_MAX == 1
							p_messages[j] = -1;
							for (int k = 0; k < nStates[s]; k++)
							{
								mesg = incoming[k] * p_edgePot[0];
								if (mesg > p_messages[j])
									p_messages[j] = mesg;
								p_edgePot += maxState;
							}
#else
							p_messages[j] = 0;
							for (int k = 0; k < nStates[s]; k++)
							{
								p_messages[j] += incoming[k] * p_edgePot[0];
								p_edgePot += maxState;
							}
#endif
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

/* Loopy BP */

#if CRF_MAX == 1
void CRF::LoopyBP_max(double *messages_1, double *messages_2, int maxIter, double cutoff, int verbose)
#else
void CRF::LoopyBP(double *messages_1, double *messages_2, int maxIter, double cutoff, int verbose)
#endif
{
	double *old_messages_1 = (double *) R_alloc(maxState * nEdges, sizeof(double));
	double *old_messages_2 = (double *) R_alloc(maxState * nEdges, sizeof(double));
	for (int i = 0; i < maxState * nEdges; i++)
		messages_1[i] = messages_2[i] = old_messages_1[i] = old_messages_2[i] = 0;

	double *incoming = (double *) R_alloc(maxState, sizeof(double));

	int s, r, e, n;
	double sumMesg, *p_nodePot, *p_edgePot, *p0_edgePot, *p_messages;
#if CRF_MAX == 1
	double mesg;
#endif

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
#if CRF_MAX == 1
						p_messages[j] = -1;
						for (int k = 0; k < nStates[s]; k++)
						{
							mesg = incoming[k] * p_edgePot[k];
							if (mesg > p_messages[j])
								p_messages[j] = mesg;
						}
#else
						p_messages[j] = 0;
						for (int k = 0; k < nStates[s]; k++)
							p_messages[j] += incoming[k] * p_edgePot[k];
#endif
						sumMesg += p_messages[j];
					}
				}
				else
				{
					p_messages = messages_1 + maxState * e;
					for (int j = 0; j < nStates[r]; j++)
					{
						p_edgePot = p0_edgePot++;
#if CRF_MAX == 1
						p_messages[j] = -1;
						for (int k = 0; k < nStates[s]; k++)
						{
							mesg = incoming[k] * p_edgePot[0];
							if (mesg > p_messages[j])
								p_messages[j] = mesg;
							p_edgePot += maxState;
						}
#else
						p_messages[j] = 0;
						for (int k = 0; k < nStates[s]; k++)
						{
							p_messages[j] += incoming[k] * p_edgePot[0];
							p_edgePot += maxState;
						}
#endif
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

#endif // CRF_MAX
