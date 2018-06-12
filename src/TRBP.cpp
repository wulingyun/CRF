#include "CRF.h"
#include <Rmath.h>

/* Tree-Reweighted BP */

void CRF::TRBP(double *mu, double **scaleEdgePot, int maxIter, double cutoff, int verbose, bool maximize)
{
	swap2(edgePot, scaleEdgePot);

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
	double *msg, sumBel;

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
			/* gather incoming messages */

			sumBel = 0;
			for (int i = 0; i < nStates[s]; i++)
			  sumBel += NodeBel(s, i) = NodePot(s, i);
			for (int i = 0; i < nStates[s]; i++)
			  NodeBel(s, i) /= sumBel;
			for (int i = 0; i < nAdj[s]; i++)
			{
				e = AdjEdges(s, i);
				if (EdgesBegin(e) == s)
					msg = old_messages[0][e];
				else
					msg = old_messages[1][e];
				sumBel = 0;
				for (int k = 0; k < nStates[s]; k++)
				  sumBel += NodeBel(s, k) *= R_pow(msg[k], mu[e]);
				for (int k = 0; k < nStates[s]; k++)
				  NodeBel(s, k) /= sumBel;
			}

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
			Rprintf("TRBP: Iteration %d, Difference = %f\n", iter, difference);
		if (difference <= cutoff)
			break;
	}

	if (difference > cutoff)
		warning("Tree-Reweighted BP did not converge in %d iterations! (diff = %f)", maxIter, difference);

	swap2(edgePot, scaleEdgePot);
}

/* Edge beliefs */

void CRF::TRBP_Messages2EdgeBel(double *mu, double **scaleEdgePot)
{
	for (int i=0; i < nEdges; i++)
	{
		for (int j = 0; j < nEdgeStates[i]; j++)
			edgeBel[i][j] = scaleEdgePot[i][j];
	}

	int n1, n2;
	double bel, sumBel;
	for (int i = 0; i < nEdges; i++)
	{
		n1 = EdgesBegin(i);
		n2 = EdgesEnd(i);
		for (int j = 0; j < nStates[n1]; j++)
		{
			bel = messages[0][i][j] == 0 ? 0 : NodeBel(n1, j) / messages[0][i][j];
			for (int k = 0; k < nStates[n2]; k++)
				EdgeBel(i, j, k) *= bel;
		}
		for (int j = 0; j < nStates[n2]; j++)
		{
			bel = messages[1][i][j] == 0 ? 0 : NodeBel(n2, j) / messages[1][i][j];
			for (int k = 0; k < nStates[n1]; k++)
				EdgeBel(i, k, j) *= bel;
		}

		sumBel = 0;
		for (int j = 0; j < nStates[n2]; j++)
		{
			for (int k = 0; k < nStates[n1]; k++)
				sumBel += EdgeBel(i, k, j);
		}
		for (int j = 0; j < nStates[n2]; j++)
		{
			for (int k = 0; k < nStates[n1]; k++)
				EdgeBel(i, k, j) /= sumBel;
		}
	}
}

/* Bethe free energy */

void CRF::TRBP_BetheFreeEnergy(double *mu)
{
	double nodeEnergy, nodeEntropy, edgeEnergy, edgeEntropy;
	nodeEnergy = nodeEntropy = edgeEnergy = edgeEntropy = 0;

	double entropy, bel, sum_mu;
	for (int i = 0; i < nNodes; i++)
	{
		entropy = 0;
		for (int j = 0; j < nStates[i]; j++)
		{
			bel = NodeBel(i, j);
			if (bel > 0)
			{
				nodeEnergy -= bel * log(NodePot(i, j));
				entropy += bel * log(bel);
			}
		}
		sum_mu = 0;
		for (int j = 0; j < nAdj[i]; j++)
			sum_mu += mu[AdjEdges(i, j)];
		nodeEntropy += (sum_mu - 1) * entropy;
	}

	int n1, n2;
	for (int i = 0; i < nEdges; i++)
	{
		entropy = 0;
		n1 = EdgesBegin(i);
		n2 = EdgesEnd(i);
		for (int j = 0; j < nStates[n2]; j++)
		{
			for (int k = 0; k < nStates[n1]; k++)
			{
				bel = EdgeBel(i, k, j);
				if (bel > 0)
				{
					edgeEnergy -= bel * log(EdgePot(i, k, j));
					entropy -= bel * log(bel);
				}
			}
		}
		edgeEntropy += mu[i] * entropy;
	}

	*logZ = - nodeEnergy + nodeEntropy - edgeEnergy + edgeEntropy;
}

/* initialize TRBP parameters */

void CRF::TRBP_Init(double *mu, double **scaleEdgePot)
{
	/* Calculate Tree Weights */

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

	/* scale edge potentials */
	
	double inv_mu;
	for (int i = 0; i < nEdges; i++)
	{
		inv_mu = 1/mu[i];
		for (int j = 0; j < nEdgeStates[i]; j++)
			scaleEdgePot[i][j] = R_pow(edgePot[i][j], inv_mu);
	}
}
