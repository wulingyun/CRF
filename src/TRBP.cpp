#include "CRF.h"
#include <Rmath.h>

/* Tree-Reweighted BP */

void CRF::TRBP(double *messages_1, double *messages_2, double *mu, double *scaleEdgePot, int maxIter, double cutoff, int verbose, bool maximize)
{
	double *old_messages_1 = (double *) R_alloc(maxState * nEdges, sizeof(double));
	double *old_messages_2 = (double *) R_alloc(maxState * nEdges, sizeof(double));
	for (int i = 0; i < maxState * nEdges; i++)
		old_messages_1[i] = old_messages_2[i] = 0;

	double *incoming = (double *) R_alloc(maxState, sizeof(double));
	double *outgoing = (double *) R_alloc(maxState, sizeof(double));

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

			p_nodePot = nodePot + s;
			for (int i = 0; i < nStates[s]; i++)
			{
				incoming[i] = p_nodePot[0];
				p_nodePot += nNodes;
			}
			for (int i = 0; i < nAdj[s]; i++)
			{
				e = adjEdges[s][i] - 1;
				if (edges[e] - 1 == s)
					p_messages = old_messages_1;
				else
					p_messages = old_messages_2;
				p_messages += maxState * e;
				for (int k = 0; k < nStates[s]; k++)
					incoming[k] *= R_pow(p_messages[k], mu[e]);
			}

			/* send messages */

			for (int i = 0; i < nAdj[s]; i++)
			{
				r = adjNodes[s][i] - 1;
				e = adjEdges[s][i] - 1;

				if (edges[e] - 1 == s)
					p_messages = old_messages_1;
				else
					p_messages = old_messages_2;
				p_messages += maxState * e;
				for (int k = 0; k < nStates[s]; k++)
					outgoing[k] = p_messages[k] == 0 ? 0 : incoming[k] / p_messages[k];

				sumMesg = 0;
				p0_edgePot = scaleEdgePot + maxState * maxState * e;
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
								mesg = outgoing[k] * p_edgePot[k];
								if (mesg > p_messages[j])
									p_messages[j] = mesg;
							}
						}
						else
						{
							for (int k = 0; k < nStates[s]; k++)
								p_messages[j] += outgoing[k] * p_edgePot[k];
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
								mesg = outgoing[k] * p_edgePot[0];
								if (mesg > p_messages[j])
									p_messages[j] = mesg;
								p_edgePot += maxState;
							}
						}
						else
						{
							for (int k = 0; k < nStates[s]; k++)
							{
								p_messages[j] += outgoing[k] * p_edgePot[0];
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
	for (int i = 0; i < nNodes; i++)
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
			if (mu[i] <= 0)
			{
				loop = 1;
				break;
			}
	}
	PutRNGstate();

	for (int i = 0; i < nEdges; i++)
		mu[i] /= n;
}

/* Node beliefs */

void CRF::TRBP_Message2NodeBelief(double *messages_1, double *messages_2, double *mu)
{
	for (int i = 0; i < length(_nodePot); i++)
		nodeBel[i] = nodePot[i];

	int n1, n2;
	double sumBel, *p_nodeBel;
	double *p1_messages = messages_1;
	double *p2_messages = messages_2;
	for (int i = 0; i < nEdges; i++)
	{
		n1 = edges[i] - 1;
		n2 = edges[i + nEdges] - 1;
		p_nodeBel = nodeBel + n1;
		for (int j = 0; j < nStates[n1]; j++)
		{
			p_nodeBel[0] *= R_pow(p1_messages[j], mu[i]);
			p_nodeBel += nNodes;
		}
		p_nodeBel = nodeBel + n2;
		for (int j = 0; j < nStates[n2]; j++)
		{
			p_nodeBel[0] *= R_pow(p2_messages[j], mu[i]);
			p_nodeBel += nNodes;
		}
		p1_messages += maxState;
		p2_messages += maxState;
	}

	for (int i = 0; i < nNodes; i++)
	{
		sumBel = 0;
		p_nodeBel = nodeBel + i;
		for (int j = 0; j < nStates[i]; j++)
		{
			sumBel += p_nodeBel[0];
			p_nodeBel += nNodes;
		}
		p_nodeBel = nodeBel + i;
		for (int j = 0; j < nStates[i]; j++)
		{
			p_nodeBel[0] /= sumBel;
			p_nodeBel += nNodes;
		}
	}
}

/* Edge beliefs */

void CRF::TRBP_Message2EdgeBelief(double *messages_1, double *messages_2, double *mu, double *scaleEdgePot)
{
	for (int i = 0; i < length(_edgePot); i++)
		edgeBel[i] = scaleEdgePot[i];

	int n1, n2;
	double bel, sumBel, *p_nodeBel, *p_edgeBel;
	double *p0_edgeBel = edgeBel;
	double *p1_messages = messages_1;
	double *p2_messages = messages_2;
	for (int i = 0; i < nEdges; i++)
	{
		n1 = edges[i] - 1;
		n2 = edges[i + nEdges] - 1;
		p_nodeBel = nodeBel + n1;
		for (int j = 0; j < nStates[n1]; j++)
		{
			bel = p1_messages[j] == 0 ? 0 : p_nodeBel[0] / p1_messages[j];
			p_edgeBel = p0_edgeBel + j;
			for (int k = 0; k < nStates[n2]; k++)
			{
				p_edgeBel[0] *= bel;
				p_edgeBel += maxState;
			}
			p_nodeBel += nNodes;
		}
		p_nodeBel = nodeBel + n2;
		p_edgeBel = p0_edgeBel;
		for (int j = 0; j < nStates[n2]; j++)
		{
			bel = p2_messages[j] == 0 ? 0 : p_nodeBel[0] / p2_messages[j];
			for (int k = 0; k < nStates[n1]; k++)
				p_edgeBel[k] *= bel;
			p_nodeBel += nNodes;
			p_edgeBel += maxState;
		}

		sumBel = 0;
		p_edgeBel = p0_edgeBel;
		for (int j = 0; j < nStates[n2]; j++)
		{
			for (int k = 0; k < nStates[n1]; k++)
				sumBel += p_edgeBel[k];
			p_edgeBel += maxState;
		}
		p_edgeBel = p0_edgeBel;
		for (int j = 0; j < nStates[n2]; j++)
		{
			for (int k = 0; k < nStates[n1]; k++)
				p_edgeBel[k] /= sumBel;
			p_edgeBel += maxState;
		}

		p0_edgeBel += maxState * maxState;
		p1_messages += maxState;
		p2_messages += maxState;
	}
}

/* Bethe free energy */

void CRF::TRBP_BetheFreeEnergy(double *mu)
{
	double nodeEnergy, nodeEntropy, edgeEnergy, edgeEntropy;
	nodeEnergy = nodeEntropy = edgeEnergy = edgeEntropy = 0;

	double entropy, sum_mu;
	double *p_nodeBel, *p_nodePot;
	for (int i = 0; i < nNodes; i++)
	{
		entropy = 0;
		p_nodeBel = nodeBel + i;
		p_nodePot = nodePot + i;
		for (int j = 0; j < nStates[i]; j++)
		{
			if (p_nodeBel[0] > 0)
			{
				nodeEnergy -= p_nodeBel[0] * log(p_nodePot[0]);
				entropy += p_nodeBel[0] * log(p_nodeBel[0]);
			}
			p_nodeBel += nNodes;
			p_nodePot += nNodes;
		}
		sum_mu = 0;
		for (int j = 0; j < nAdj[i]; j++)
			sum_mu += mu[adjEdges[i][j] - 1];
		nodeEntropy += (sum_mu - 1) * entropy;
	}

	int n1, n2;
	double *p_edgeBel, *p0_edgeBel = edgeBel;
	double *p_edgePot, *p0_edgePot = edgePot;
	for (int i = 0; i < nEdges; i++)
	{
		entropy = 0;
		n1 = edges[i] - 1;
		n2 = edges[i + nEdges] - 1;
		p_edgeBel = p0_edgeBel;
		p_edgePot = p0_edgePot;
		for (int j = 0; j < nStates[n2]; j++)
		{
			for (int k = 0; k < nStates[n1]; k++)
			{
				if (p_edgeBel[k] > 0)
				{
					edgeEnergy -= p_edgeBel[k] * log(p_edgePot[k]);
					entropy -= p_edgeBel[k] * log(p_edgeBel[k]);
				}
			}
			p_edgeBel += maxState;
			p_edgePot += maxState;
		}
		p0_edgeBel += maxState * maxState;
		p0_edgePot += maxState * maxState;
		edgeEntropy += mu[i] * entropy;
	}

	*logZ = - nodeEnergy + nodeEntropy - edgeEnergy + edgeEntropy;
}

void CRF::TRBP_ScaleEdgePot(double *mu, double *scaleEdgePot)
{
	double inv_mu;
	int ss = maxState * maxState, n = 0;
	for (int i = 0; i < nEdges; i++)
	{
		inv_mu = 1/mu[i];
		for (int j = 0; j < ss; j++)
		{
			scaleEdgePot[n] = R_pow(edgePot[n], inv_mu);
			n++;
		}
	}
}
