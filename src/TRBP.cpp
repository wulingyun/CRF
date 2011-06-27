#include "CRF.h"
#include <Rmath.h>

/* Tree-Reweighted BP */

void CRF::TRBP(double *messages_1, double *messages_2, double *mu, int maxIter, double cutoff, int verbose, bool maximize)
{
	double *old_messages_1 = (double *) R_alloc(maxState * nEdges, sizeof(double));
	double *old_messages_2 = (double *) R_alloc(maxState * nEdges, sizeof(double));
	for (int i = 0; i < maxState * nEdges; i++)
		messages_1[i] = messages_2[i] = old_messages_1[i] = old_messages_2[i] = 0;

	double *inv_mu = (double *) R_alloc(nEdges, sizeof(double));
	for (int i = 0; i < nEdges; i++)
		inv_mu[i] = 1/mu[i];

	double *incoming = (double *) R_alloc(maxState, sizeof(double));

	int s, r, e, n;
	double mesg, sumMesg, *p_nodePot, *p_edgePot, *p0_edgePot, *p_messages;

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
					e = adjEdges[s][j] - 1;
					if (edges[e] - 1 == s)
						p_messages = old_messages_1;
					else
						p_messages = old_messages_2;
					p_messages += maxState * e;
					if (j != i)
						for (int k = 0; k < nStates[s]; k++)
							incoming[k] *= R_pow(p_messages[k], mu[e]);
					else
						for (int k = 0; k < nStates[s]; k++)
							incoming[k] /= R_pow(p_messages[k], 1-mu[e]);
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
						p_messages[j] = 0;
						if (maximize)
						{
							for (int k = 0; k < nStates[s]; k++)
							{
								mesg = incoming[k] * R_pow(p_edgePot[k], inv_mu[e]);
								if (mesg > p_messages[j])
									p_messages[j] = mesg;
							}
						}
						else
						{
							for (int k = 0; k < nStates[s]; k++)
								p_messages[j] += incoming[k] * R_pow(p_edgePot[k], inv_mu[e]);
						}
						sumMesg += p_messages[j];
					}
				}
				else
				{
					p_messages = messages_1 + maxState * e;
					for (int j = 0; j < nStates[r]; j++)
					{
						p_edgePot = p0_edgePot++;
						p_messages[j] = 0;
						if (maximize)
						{
							for (int k = 0; k < nStates[s]; k++)
							{
								mesg = incoming[k] * R_pow(p_edgePot[0], inv_mu[e]);
								if (mesg > p_messages[j])
									p_messages[j] = mesg;
								p_edgePot += maxState;
							}
						}
						else
						{
							for (int k = 0; k < nStates[s]; k++)
							{
								p_messages[j] += incoming[k] * R_pow(p_edgePot[0], inv_mu[e]);
								p_edgePot += maxState;
							}
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
			Rprintf("TRBP: Iteration %d, Difference = %f\n", iter, difference);
		if (difference <= cutoff)
			break;
	}

	if (difference > cutoff)
		warning("Tree-Reweighted BP did not converge in %d iterations! (diff = %f)", maxIter, difference);
}

/* Minimum Weight Spanning Tree using Kruskal algorithm */

void MinSpanTree(int *tree, int nNodes, int nEdges, int *edges, double *costs)
{
	int *index = (int *) R_alloc(nEdges, sizeof(int));
	for (int i = 0; i < nEdges; i++)
	{
		tree[i] = 0;
		index[i] = i;
	}
	rsort_with_index(costs, index, nEdges);

	int *label = (int *) R_alloc(nNodes, sizeof(int));
	for (int i = 0; i < nEdges; i++)
		label[i] = i;

	int n = 0, n1, n2;
	for (int i = 0; i < nEdges; i++)
	{
		n1 = edges[index[i]] - 1;
		n2 = edges[index[i] + nEdges] - 1;
		if (label[n1] != label[n2])
		{
			for (int j = 0; j < nEdges; j++)
				if (label[j] == label[n2])
					label[j] = label[n1];
			tree[index[i]] = 1;
			if (++n >= nNodes - 1)
				break;
		}
	}
}

/* Calculate Tree Weights */

void CRF::TRBP_Weights(double *mu)
{
	for (int i = 0; i < nEdges; i++)
		mu[i] = 0;

	int *tree = (int *) R_alloc(nEdges, sizeof(int));
	double *costs = (double *) R_alloc(nEdges, sizeof(double));
	int n = 0, loop = 1;

	GetRNGstate();
	while (loop)
	{
		for (int i = 0; i < nEdges; i++)
			costs[i] = unif_rand();

		MinSpanTree(tree, nNodes, nEdges, edges, costs);

		for (int i = 0; i < nEdges; i++)
			if (tree[i])
				mu[i]++;
		n++;

		loop = 0;
		for (int i = 0; i < nEdges; i++)
			if (mu[i] == 0)
			{
				loop = 1;
				break;
			}
	}
	PutRNGstate();

	for (int i = 0; i < nEdges; i++)
		mu[i] /= n;
}
