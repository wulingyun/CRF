#include "CRF.h"

/* BP messages init */

void CRF::MessagesInit()
{
	int dim[] = {2, nEdges, maxState};
	messages = (double ***) allocArray<double, 3>(dim);
}

/* Node beliefs */

void CRF::Messages2NodeBel()
{
	for (int i = 0; i < length(_nodePot); i++)
		nodeBel[i] = nodePot[i];

	int n1, n2;
	double sumBel;
	for (int i = 0; i < nEdges; i++)
	{
		n1 = EdgesBegin(i);
		n2 = EdgesEnd(i);
		for (int j = 0; j < nStates[n1]; j++)
			NodeBel(n1, j) *= messages[0][i][j];
		for (int j = 0; j < nStates[n2]; j++)
			NodeBel(n2, j) *= messages[1][i][j];
	}

	for (int i = 0; i < nNodes; i++)
	{
		sumBel = 0;
		for (int j = 0; j < nStates[i]; j++)
			sumBel += NodeBel(i, j);
		for (int j = 0; j < nStates[i]; j++)
			NodeBel(i, j) /= sumBel;
	}
}

/* Edge beliefs */

void CRF::Messages2EdgeBel()
{
	for (int i = 0; i < nEdges; i++)
	{
		for (int j = 0; j < nEdgeStates[i]; j++)
			edgeBel[i][j] = edgePot[i][j];
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

/* Decoding by max of marginals */

void CRF::MaxOfMarginals()
{
	double maxBel;
	for (int i = 0; i < nNodes; i++)
	{
		maxBel = -1;
		for (int j = 0; j < nStates[i]; j++)
		{
			if (NodeBel(i, j) > maxBel)
			{
				maxBel = NodeBel(i, j);
				labels[i] = j;
			}
		}
	}

	for (int i = 0; i < nNodes; i++)
		labels[i]++;
}

/* Bethe free energy */

void CRF::BetheFreeEnergy()
{
	double nodeEnergy, nodeEntropy, edgeEnergy, edgeEntropy;
	nodeEnergy = nodeEntropy = edgeEnergy = edgeEntropy = 0;

	double entropy, bel;
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
		nodeEntropy += (nAdj[i] - 1) * entropy;
	}

	int n1, n2;
	for (int i = 0; i < nEdges; i++)
	{
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
					edgeEntropy -= bel * log(bel);
				}
			}
		}
	}

	*logZ = - nodeEnergy + nodeEntropy - edgeEnergy + edgeEntropy;
}
