#include "CRF.h"

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

	*logZ = - (nodeEnergy + edgeEnergy - nodeEntropy - edgeEntropy);
}

/* compute sum-product messages */

double *CRF::ComputeMessagesSum(int s, int r, int e, double *outgoing, double ***old_messages, double ***new_messages)
{
	double *msg, sumMsg = 0;
	if (EdgesBegin(e) == s)
	{
		for (int j = 0; j < nStates[s]; j++)
			outgoing[j] = old_messages[0][e][j] == 0 ? 0 : NodeBel(s, j) / old_messages[0][e][j];
		msg = new_messages[1][e];
		for (int j = 0; j < nStates[r]; j++)
		{
			msg[j] = 0;
			for (int k = 0; k < nStates[s]; k++)
				msg[j] += outgoing[k] * EdgePot(e, k, j);
			sumMsg += msg[j];
		}
	}
	else
	{
		for (int j = 0; j < nStates[s]; j++)
			outgoing[j] = old_messages[1][e][j] == 0 ? 0 : NodeBel(s, j) / old_messages[1][e][j];
		msg = new_messages[0][e];
		for (int j = 0; j < nStates[r]; j++)
		{
			msg[j] = 0;
			for (int k = 0; k < nStates[s]; k++)
				msg[j] += outgoing[k] * EdgePot(e, j, k);
			sumMsg += msg[j];
		}
	}
	for (int j = 0; j < nStates[r]; j++)
		msg[j] /= sumMsg;
	return msg;
}

/* compute max-product messages */

double *CRF::ComputeMessagesMax(int s, int r, int e, double *outgoing, double ***old_messages, double ***new_messages)
{
	double *msg, m, sumMsg = 0;
	if (EdgesBegin(e) == s)
	{
		for (int j = 0; j < nStates[s]; j++)
			outgoing[j] = old_messages[0][e][j] == 0 ? 0 : NodeBel(s, j) / old_messages[0][e][j];
		msg = new_messages[1][e];
		for (int j = 0; j < nStates[r]; j++)
		{
			msg[j] = 0;
			for (int k = 0; k < nStates[s]; k++)
			{
				m = outgoing[k] * EdgePot(e, k, j);
				if (m > msg[j])
					msg[j] = m;
			}
			sumMsg += msg[j];
		}
	}
	else
	{
		for (int j = 0; j < nStates[s]; j++)
			outgoing[j] = old_messages[1][e][j] == 0 ? 0 : NodeBel(s, j) / old_messages[1][e][j];
		msg = new_messages[0][e];
		for (int j = 0; j < nStates[r]; j++)
		{
			msg[j] = 0;
			for (int k = 0; k < nStates[s]; k++)
			{
				m = outgoing[k] * EdgePot(e, j, k);
				if (m > msg[j])
					msg[j] = m;
			}
			sumMsg += msg[j];
		}
	}
	for (int j = 0; j < nStates[r]; j++)
		msg[j] /= sumMsg;
	return msg;
}

void CRF::GatherIncomingMessages(int s, double ***old_messages)
{
  double sumBel, *msg;
  int e;
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
      sumBel += NodeBel(s, k) *= msg[k];
    for (int k = 0; k < nStates[s]; k++)
      NodeBel(s, k) /= sumBel;
  }
}
