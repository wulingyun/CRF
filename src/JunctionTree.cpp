#include "CRF.h"

/* build junction tree */

JunctionTree::JunctionTree(CRF &crf)
: original(crf)
{
	nNodes = original.nNodes;
	nEdges = original.nEdges;
	nStates = original.nStates;

	nClusters = 0;
	nClusterNodes = (int *) R_allocVector<int>(nNodes);

	/* triangulation using minimum fill-in heuristic */

	int **cliqueNodes = (int **) C_allocArray<int>(nNodes, nNodes);
	int **adj = (int **) C_allocArray<int>(nNodes, nNodes);
	int **neighbors = (int **) C_allocArray<int>(nNodes, nNodes);
	int *nNeighbors = (int *) C_allocVector<int>(nNodes);
	int *nMissingEdges = (int *) C_allocVector<int>(nNodes);
	int *overlap = (int *) C_allocVector<int>(nNodes);

	for (int i = 0; i < nNodes; i++)
	{
		for (int j = 0; j < nNodes; j++)
			adj[i][j] = 0;
		nNeighbors[i] = 0;
		nClusterNodes[i] = 0;
	}

	int n1, n2;
	for (int i = 0; i < nEdges; i++)
	{
		n1 = original.EdgesBegin(i);
		n2 = original.EdgesEnd(i);
		Insert(neighbors[n1], nNeighbors[n1], n2);
		Insert(neighbors[n2], nNeighbors[n2], n1);
		adj[n1][n2] = adj[n2][n1] = 1;
	}

	for (int i = 0; i < nNodes; i++)
	{
		nMissingEdges[i] = 0;
		for (int k1 = 0; k1 < nNeighbors[i]-1; k1++)
		{
			n1 = neighbors[i][k1];
			for (int k2 = k1+1; k2 < nNeighbors[i]; k2++)
			{
				n2 = neighbors[i][k2];
				if (adj[n1][n2] == 0)
					nMissingEdges[i]++;
			}
		}
	}

	int n, m, maxMissingEdges = nNodes * nNodes;
	int treeWidth = 0;
	while (1)
	{
		n = -1;
		m = maxMissingEdges;
		for (int i = 0; i < nNodes; i++)
		{
			if (nMissingEdges[i] >= 0 && nMissingEdges[i] < m)
			{
				n = i;
				m = nMissingEdges[i];
			}
		}
		if (n < 0)
			break;

		for (int j1 = 0; j1 < nNeighbors[n]-1; j1++)
		{
			n1 = neighbors[n][j1];
			for (int j2 = j1+1; j2 < nNeighbors[n]; j2++)
			{
				n2 = neighbors[n][j2];
				if (adj[n1][n2] == 0)
				{
					m = Intersection(overlap, neighbors[n1], nNeighbors[n1], neighbors[n2], nNeighbors[n2]);
					for (int k = 0; k < m; k++)
						nMissingEdges[overlap[k]]--;
					nMissingEdges[n1] += nNeighbors[n1]-m;
					nMissingEdges[n2] += nNeighbors[n2]-m;
					Insert(neighbors[n1], nNeighbors[n1], n2);
					Insert(neighbors[n2], nNeighbors[n2], n1);
					adj[n1][n2] = adj[n2][n1] = 1;
				}
			}
		}
		for (int j1 = 0; j1 < nNeighbors[n]; j1++)
		{
			n1 = neighbors[n][j1];
			m = Intersection(overlap, neighbors[n1], nNeighbors[n1], neighbors[n], nNeighbors[n]);
			nMissingEdges[n1] -= nNeighbors[n1]-1-m;
			Remove(neighbors[n1], nNeighbors[n1], n);
			Insert(cliqueNodes[nClusters], nClusterNodes[nClusters], n1);
		}
		nMissingEdges[n] = -1;
		nNeighbors[n] = 0;
		Insert(cliqueNodes[nClusters], nClusterNodes[nClusters], n);
		if (treeWidth < nClusterNodes[nClusters])
			treeWidth = nClusterNodes[nClusters];
		else if (treeWidth > nClusterNodes[nClusters])
		{
			for (int i = 0; i < nClusters; i++)
			{
				if (nClusterNodes[i] > nClusterNodes[nClusters] && cliqueNodes[i][0] <= cliqueNodes[nClusters][0] &&
					cliqueNodes[i][nClusterNodes[i]-1] >= cliqueNodes[nClusters][nClusterNodes[nClusters]-1])
				{
					m = Intersection(overlap, cliqueNodes[i], nClusterNodes[i], cliqueNodes[nClusters], nClusterNodes[nClusters]);
					if (m == nClusterNodes[nClusters])
					{
						nClusterNodes[nClusters] = 0;
						break;
					}
				}
			}
		}
		if (nClusterNodes[nClusters] > 0)
			nClusters++;
	}
	C_freeArray<int, 2>(adj);
	C_freeArray<int, 2>(neighbors);
	C_freeVector(nNeighbors);
	C_freeVector(nMissingEdges);

	/* cluster information */

	int **cliqueEdges = C_allocArray<int>(nClusters, nEdges);
	nClusterEdges = (int *) R_allocVector<int>(nClusters);
	for (int i = 0; i < nClusters; i++)
	{
		for (int j = 0; j < nNodes; j++)
			overlap[j] = 0;
		for (int j = 0; j < nClusterNodes[i]; j++)
			overlap[cliqueNodes[i][j]] = 1;

		nClusterEdges[i] = 0;
		for (int j = 0; j < nClusterNodes[i]; j++)
		{
			n = cliqueNodes[i][j];
			for (int k = 0; k < original.nAdj[n]; k++)
			{
				m = original.AdjEdges(n, k);
				if (original.AdjNodes(n, k) > n && overlap[original.EdgesEnd(m)])
					cliqueEdges[i][nClusterEdges[i]++] = m;
			}
		}
	}

	clusterNodes = (int **) R_allocArray2<int>(nClusters, nClusterNodes);
	clusterEdges = (int **) R_allocArray2<int>(nClusters, nClusterEdges);
	for (int i = 0; i < nClusters; i++)
	{
		for (int j = 0; j < nClusterNodes[i]; j++)
			clusterNodes[i][j] = cliqueNodes[i][j];
		for (int j = 0; j < nClusterEdges[i]; j++)
			clusterEdges[i][j] = cliqueEdges[i][j];
	}
	C_freeArray<int, 2>(cliqueNodes);
	C_freeArray<int, 2>(cliqueEdges);

	/* construct junction tree by max weighted spanning tree */

	int nCliqueEdges = nClusters * (nClusters - 1) / 2;
	int *tree = (int *) C_allocVector<int>(nCliqueEdges);
	int *edges = (int *) C_allocVector<int>(nCliqueEdges * 2);
	int *weights = (int *) C_allocVector<int>(nCliqueEdges);
	double *costs = (double *) C_allocVector<double>(nCliqueEdges);
	n = 0;
	for (int i = 0; i < nClusters-1; i++)
	{
		for (int j = i+1; j < nClusters; j++)
		{
			edges[n] = i;
			edges[n + nCliqueEdges] = j;
			weights[n] = Intersection(overlap, clusterNodes[i], nClusterNodes[i], clusterNodes[j], nClusterNodes[j]);
			costs[n] = -weights[n];
			n++;
		}
	}
	nSeperators = MinSpanTree(tree, nClusters, nCliqueEdges, edges, costs, 0);

	nAdj = (int *) R_allocVector<int>(nClusters);
	for (int i = 0; i < nClusters; i++)
		nAdj[i] = 0;

	seperators = (int **) R_allocArray<int>(nSeperators, 2);
	nSeperatorNodes = (int *) R_allocVector<int>(nSeperators);
	n = 0;
	for (int i = 0; i < nCliqueEdges; i++)
	{
		if (tree[i])
		{
			n1 = edges[i];
			n2 = edges[i + nCliqueEdges];
			nAdj[n1]++;
			nAdj[n2]++;
			seperators[n][0] = n1;
			seperators[n][1] = n2;
			nSeperatorNodes[n] = weights[i];
			n++;
		}
	}

	C_freeVector(overlap);
	C_freeVector(tree);
	C_freeVector(edges);
	C_freeVector(weights);
	C_freeVector(costs);

	/* auxiliary variabls for accessing junction tree */

	adjClusters = (int **) R_allocArray2<int>(nClusters, nAdj);
	adjSeperators = (int **) R_allocArray2<int>(nClusters, nAdj);
	for (int i = 0; i < nClusters; i++)
		nAdj[i] = 0;

	/* seperators (edges in junction tree) information */

	seperatorNodes = (int **) R_allocArray2<int>(nSeperators, nSeperatorNodes);
	for (int i = 0; i < nSeperators; i++)
	{
		n1 = seperators[i][0];
		n2 = seperators[i][1];
		adjClusters[n1][nAdj[n1]] = n2;
		adjClusters[n2][nAdj[n2]] = n1;
		adjSeperators[n1][nAdj[n1]] = i;
		adjSeperators[n2][nAdj[n2]] = i;
		nAdj[n1]++;
		nAdj[n2]++;
		Intersection(seperatorNodes[i], clusterNodes[n1], nClusterNodes[n1], clusterNodes[n2], nClusterNodes[n2]);
	}

	/* number of cluster states */

	nClusterStates = (int *) R_allocVector<int>(nClusters);
	for (int i = 0; i < nClusters; i++)
	{
		nClusterStates[i] = 1;
		for (int j = 0; j < nClusterNodes[i]; j++)
			nClusterStates[i] *= nStates[clusterNodes[i][j]];
	}

	/* number of seperator states */

	nSeperatorStates = (int *) R_allocVector<int>(nSeperators);
	for (int i = 0; i < nSeperators; i++)
	{
		nSeperatorStates[i] = 1;
		for (int j = 0; j < nSeperatorNodes[i]; j++)
			nSeperatorStates[i] *= nStates[seperatorNodes[i][j]];
	}

	/* beliefs of clusters and seperators */

	clusterBel = (double **) R_allocArray2<double>(nClusters, nClusterStates);
	seperatorBel = (double **) R_allocArray2<double>(nSeperators, nSeperatorStates);

	/* temporary variables for passing messages */

	masks = (int *) R_allocVector<int>(nNodes);
	states = (int *) R_allocVector<int>(nNodes);
}

int JunctionTree::States2Index(int nNodes, int *nodes, int *states)
{
	int n = nNodes-1;
	int k = states[nodes[n]];
	for (int i = n-1; i >= 0; i--)
	{
		k *= nStates[nodes[i]];
		k += states[nodes[i]];
	}
	return k;
}

int *JunctionTree::Index2States(int nNodes, int *nodes, int index, int *states)
{
	int n;
	for (int i = 0; i < nNodes-1; i++)
	{
		n = nStates[nodes[i]];
		states[nodes[i]] = index % n;
		index /= n;
	}
	states[nodes[nNodes-1]] = index;
	return states;
}

void JunctionTree::InitStateMasks(int c, int s)
{
	cid = c;
	for (int i = 0; i < nClusterNodes[cid]; i++)
		masks[clusterNodes[cid][i]] = 0;
	if (s >= 0)
	{
		sid = s;
		for (int i = 0; i < nSeperatorNodes[sid]; i++)
			masks[seperatorNodes[sid][i]] = 1;
	}
}

void JunctionTree::ResetClusterState()
{
	for (int i = 0; i < nClusterNodes[cid]; i++)
	{
		if (masks[clusterNodes[cid][i]])
			continue;
		states[clusterNodes[cid][i]] = 0;
	}
}

void JunctionTree::ResetSeperatorState()
{
	for (int i = 0; i < nSeperatorNodes[sid]; i++)
		states[seperatorNodes[sid][i]] = 0;
}

bool JunctionTree::NextClusterState()
{
	int index, n;
	for (index = 0; index < nClusterNodes[cid]; index++)
	{
		n = clusterNodes[cid][index];
		if (masks[n])
			continue;
		states[n]++;
		if (states[n] < nStates[n])
			break;
		else
			states[n] = 0;
	}
	if (index == nClusterNodes[cid])
		return 0;
	else
		return 1;
}

bool JunctionTree::NextSeperatorState()
{
	int index, n;
	for (index = 0; index < nSeperatorNodes[sid]; index++)
	{
		n = seperatorNodes[sid][index];
		states[n]++;
		if (states[n] < nStates[n])
			break;
		else
			states[n] = 0;
	}
	if (index == nSeperatorNodes[sid])
		return 0;
	else
		return 1;
}

void JunctionTree::SendMessagesFromClusterSum(int c, int s)
{
	InitStateMasks(c, s);

	double msg, sumBel = 0;
	ResetSeperatorState();
	do
	{
		msg = 0;
		ResetClusterState();
		do
		{
			msg += ClusterBel(c, states);
		}
		while (NextClusterState());

		double &bel = SeperatorBel(s, states);
		sumBel += bel = (bel == 0) ? 0 : msg / bel;
	}
	while (NextSeperatorState());

	for (int i = 0; i < nSeperatorStates[s]; i++)
		seperatorBel[s][i] /= sumBel;
}

void JunctionTree::SendMessagesFromClusterMax(int c, int s)
{
	InitStateMasks(c, s);

	double msg, sumBel = 0;
	ResetSeperatorState();
	do
	{
		msg = 0;
		ResetClusterState();
		do
		{
			double &bel = ClusterBel(c, states);
			if (bel > msg)
				msg = bel;
		}
		while (NextClusterState());

		double &bel = SeperatorBel(s, states);
		sumBel += bel = (bel == 0) ? 0 : msg / bel;
	}
	while (NextSeperatorState());

	for (int i = 0; i < nSeperatorStates[s]; i++)
		seperatorBel[s][i] /= sumBel;
}

void JunctionTree::SendMessagesFromSeperator(int s, int c)
{
	InitStateMasks(c, s);

	double msg;
	ResetSeperatorState();
	do
	{
		msg = SeperatorBel(s, states);
		ResetClusterState();
		do
		{
			ClusterBel(c, states) *= msg;
		}
		while (NextClusterState());
	}
	while (NextSeperatorState());
}

void JunctionTree::InitMessages()
{
	for (int i = 0; i < nClusters; i++)
		for (int j = 0; j < nClusterStates[i]; j++)
			clusterBel[i][j] = 1;
	for (int i = 0; i < nSeperators; i++)
		for (int j = 0; j < nSeperatorStates[i]; j++)
			seperatorBel[i][j] = 1;

	int *nodeFree = (int *) C_allocVector<int>(nNodes);
	int *edgeFree = (int *) C_allocVector<int>(nEdges);
	for (int i = 0; i < nNodes; i++)
		nodeFree[i] = 1;
	for (int i = 0; i < nEdges; i++)
		edgeFree[i] = 1;

	int n;
	for (int c = 0; c < nClusters; c++)
	{
		InitStateMasks(c);
		ResetClusterState();
		do
		{
			double &bel = ClusterBel(c, states);
			for (int i = 0; i < nClusterNodes[c]; i++)
			{
				n = clusterNodes[c][i];
				if (nodeFree[n])
					bel *= original.NodePot(n, states[n]);
			}
			for (int i = 0; i < nClusterEdges[c]; i++)
			{
				n = clusterEdges[c][i];
				if (edgeFree[n])
					bel *= original.EdgePot(n, states[original.EdgesBegin(n)], states[original.EdgesEnd(n)]);
			}
		}
		while (NextClusterState());

		for (int i = 0; i < nClusterNodes[c]; i++)
			nodeFree[clusterNodes[c][i]] = 0;
		for (int i = 0; i < nClusterEdges[c]; i++)
			edgeFree[clusterEdges[c][i]] = 0;
	}

	C_freeVector(nodeFree);
	C_freeVector(edgeFree);
}

void JunctionTree::Messages2NodeBel(bool maximize)
{
	int *nodeFree = (int *) C_allocVector<int>(nNodes);
	for (int i = 0; i < nNodes; i++)
		nodeFree[i] = 1;

	int n;
	double bel;
	for (int c = 0; c < nClusters; c++)
	{
		InitStateMasks(c);
		for (int i = 0; i < nClusterNodes[c]; i++)
		{
			n = clusterNodes[c][i];
			if (nodeFree[n])
			{
				masks[n] = 1;
				for (int j = 0; j < nStates[n]; j++)
				{
					states[n] = j;
					bel = 0;
					ResetClusterState();
					if (maximize)
					{
						do
						{
							double b = ClusterBel(c, states);
							if (b > bel)
								bel = b;
						}
						while (NextClusterState());
					}
					else
					{
						do
						{
							bel += ClusterBel(c, states);
						}
						while (NextClusterState());
					}
					original.NodeBel(n, j) = bel;
				}
				masks[n] = 0;
				nodeFree[n] = 0;
			}
		}
	}
	original.Normalize_NodeBel();

	C_freeVector(nodeFree);
}

void JunctionTree::Messages2EdgeBel()
{
	int *edgeFree = (int *) C_allocVector<int>(nEdges);
	for (int i = 0; i < nEdges; i++)
		edgeFree[i] = 1;

	int n, n1, n2;
	double bel;
	for (int c = 0; c < nClusters; c++)
	{
		InitStateMasks(c);
		for (int i = 0; i < nClusterEdges[c]; i++)
		{
			n = clusterEdges[c][i];
			if (edgeFree[n])
			{
				n1 = original.EdgesBegin(n);
				n2 = original.EdgesEnd(n);
				masks[n1] = masks[n2] = 1;
				for (int j = 0; j < nStates[n1]; j++)
				{
					states[n1] = j;
					for (int k = 0; k < nStates[n2]; k++)
					{
						states[n2] = k;
						bel = 0;
						ResetClusterState();
						do
						{
							bel += ClusterBel(c, states);
						}
						while (NextClusterState());
						original.EdgeBel(n, j, k) = bel;
					}
				}
				masks[n1] = masks[n2] = 0;
				edgeFree[n] = 0;
			}
		}
	}
	original.Normalize_EdgeBel();

	C_freeVector(edgeFree);
}

void JunctionTree::SendMessages(bool maximize)
{
	InitMessages();

	int *nWaiting = (int *) C_allocVector<int>(nClusters);
	int **waiting = (int **) C_allocArray2<int>(nClusters, nAdj);
	int *sent = (int *) C_allocVector<int>(nClusters);
	int senderHead, senderTail, nReceiver;
	int *sender = (int *) C_allocVector<int>(nClusters * 2);
	int *receiver = (int *) C_allocVector<int>(nClusters);

	senderHead = senderTail = nReceiver = 0;
	for (int i = 0; i < nClusters; i++)
	{
		nWaiting[i] = nAdj[i];
		for (int j = 0; j < nAdj[i]; j++)
			waiting[i][j] = 1;
		sent[i] = -1;
		if (nAdj[i] == 1)
			sender[senderTail++] = i;
	}

	int s, r, e, n;

	while (senderHead < senderTail)
	{
		R_CheckUserInterrupt();

		s = sender[senderHead++];
		if (sent[s] == -2) continue;

		nReceiver = 0;
		if (nWaiting[s] == 1)
		{
			for (int i = 0; i < nAdj[s]; i++)
			{
				if (waiting[s][i])
				{
					receiver[nReceiver++] = i;
					sent[s] = nAdj[s] == 1 ? -2 : i;
					break;
				}
			}
		}
		else
		{
			for (int i = 0; i < nAdj[s]; i++)
				if (sent[s] != i)
					receiver[nReceiver++] = i;
			sent[s] = -2;
		}

		/* send messages */

		for (int i = 0; i < nReceiver; i++)
		{
			n = receiver[i];
			r = adjClusters[s][n];
			e = adjSeperators[s][n];

			for (int j = 0; j < nAdj[r]; j++)
				if (adjClusters[r][j] == s)
				{
					waiting[r][j] = 0;
					nWaiting[r]--;
					break;
				}

			if (sent[r] != -2 && nWaiting[r] <= 1)
				sender[senderTail++] = r;

			if (maximize)
				SendMessagesFromClusterMax(s, e);
			else
				SendMessagesFromClusterSum(s, e);
			SendMessagesFromSeperator(e, r);
		}
	}

	C_freeVector(nWaiting);
	C_freeArray2<int>(waiting);
	C_freeVector(sent);
	C_freeVector(sender);
	C_freeVector(receiver);

	Messages2NodeBel(maximize);
}

void JunctionTree::Sample(int size)
{
	if (size <= 0)
		size = original.nSamples;
	else if (size > original.nSamples)
		original.Init_Samples(size);

	int nOrder = 0, nQueue = 0;
	int *ordered = (int *) R_allocVector<int>(nClusters);
	int *order = (int *) R_allocVector<int>(nClusters);
	int *root = (int *) R_allocVector<int>(nClusters);
	int *queue = (int *) R_allocVector<int>(nClusters);

	for (int i = 0; i < nClusters; i++)
		ordered[i] = root[i] = 0;

	int n1, n2;
	for (int i = 0; i < nClusters; i++)
	{
		if (ordered[i])
			continue;

		ordered[i] = root[i] = 1;
		order[nOrder++] = i;

		queue[nQueue++] = i;
		while (nQueue > 0)
		{
			n1 = queue[--nQueue];
			for (int j = 0; j < nAdj[n1]; j++)
			{
				n2 = adjClusters[n1][j];
				if (ordered[n2])
					continue;

				ordered[n2] = 1;
				order[nOrder++] = n2;

				queue[nQueue++] = n2;
			}
		}
	}

	int maxState = 0;
	for (int i = 0; i < nClusters; i++)
		if (maxState < nClusterStates[i])
			maxState = nClusterStates[i];

	int *freeNodes = (int *) R_allocVector<int>(nNodes);
	double sumProb, *prob = (double *) R_allocVector<double>(maxState);

	int c, n, index, nFreeNodes, nProb;
	GetRNGstate();
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < nNodes; j++)
			states[j] = -1;

		for (int j = 0; j < nClusters; j++)
		{
			c = order[j];
			if (root[j])
			{
				nProb = nClusterStates[c];
				sumProb = 0;
				for (int k = 0; k < nProb; k++)
					sumProb += prob[k] = clusterBel[c][k];

				for (int k = 0; k < nProb; k++)
					prob[k] /= sumProb;
				index = SampleFrom(nProb, prob);
				Index2States(nClusterNodes[c], clusterNodes[c], index, states);
			}
			else
			{
				InitStateMasks(c);
				nFreeNodes = 0;
				for (int k = 0; k < nClusterNodes[c]; k++)
				{
					n = clusterNodes[c][k];
					if (states[n] >= 0)
						masks[n] = 1;
					else
						freeNodes[nFreeNodes++] = n;
				}

				nProb = 0;
				sumProb = 0;
				ResetClusterState();
				do
				{
					sumProb += prob[nProb++] = ClusterBel(c, states);
				}
				while (NextClusterState());

				for (int k = 0; k < nProb; k++)
					prob[k] /= sumProb;
				index = SampleFrom(nProb, prob);
				Index2States(nFreeNodes, freeNodes, index, states);
			}
		}

		for (int j = 0; j < nNodes; j++)
			original.Samples(i, j) = states[j] + 1;
	}
	PutRNGstate();
}
