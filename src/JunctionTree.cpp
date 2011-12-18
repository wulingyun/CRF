#include "CRF.h"

/* build junction tree */

JunctionTree::JunctionTree(CRF &crf)
: original(crf)
{
	int nNodes = original.nNodes;
	clusterSize = (int *) R_allocVector<int>(nNodes);
	nClusters = 0;

	int **cliques = (int **) C_allocArray<int>(nNodes, nNodes);
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
		clusterSize[i] = 0;
	}

	int n1, n2;
	for (int i = 0; i < original.nEdges; i++)
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
				{
					nMissingEdges[i]++;
				}
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
			Insert(cliques[nClusters], clusterSize[nClusters], n1);
		}
		nMissingEdges[n] = -1;
		nNeighbors[n] = 0;
		Insert(cliques[nClusters], clusterSize[nClusters], n);
		if (treeWidth < clusterSize[nClusters])
			treeWidth = clusterSize[nClusters];
		else if (treeWidth > clusterSize[nClusters])
		{
			for (int i = 0; i < nClusters; i++)
			{
				if (clusterSize[i] > clusterSize[nClusters] && cliques[i][0] <= cliques[nClusters][0] && cliques[i][clusterSize[i]-1] >= cliques[nClusters][clusterSize[nClusters]-1])
				{
					m = Intersection(overlap, cliques[i], clusterSize[i], cliques[nClusters], clusterSize[nClusters]);
					if (m == clusterSize[nClusters])
					{
						clusterSize[nClusters] = 0;
						break;
					}
				}
			}
		}
		if (clusterSize[nClusters] > 0)
			nClusters++;
	}
	treeWidth--;
	Rprintf("%d, %d\n", nClusters, treeWidth);

	clusters = (int **) R_allocArray2<int>(nClusters, clusterSize);
	for (int i = 0; i < nClusters; i++)
		for (int j = 0; j < clusterSize[i]; j++)
			clusters[i][j] = cliques[i][j];
	C_freeArray<int, 2>(cliques);
	C_freeArray<int, 2>(adj);
	C_freeArray<int, 2>(neighbors);
	C_freeVector(nNeighbors);
	C_freeVector(nMissingEdges);

	m = nClusters * (nClusters - 1);
	int *tree = (int *) C_allocVector<int>(m);
	int *edges = (int *) C_allocVector<int>(m * 2);
	int *weights = (int *) C_allocVector<int>(m);
	double *costs = (double *) C_allocVector<double>(m);
	n = 0;
	for (int i = 0; i < nClusters-1; i++)
	{
		for (int j = i+1; j < nClusters; j++)
		{
			edges[n] = i;
			edges[n + m] = j;
			if (clusters[i][0] <= clusters[nClusters][clusterSize[nClusters]-1] && clusters[i][clusterSize[i]-1] >= clusters[nClusters][0])
				weights[n] = Intersection(overlap, clusters[i], clusterSize[i], clusters[j], clusterSize[j]);
			else
				weights[n] = 0;
			costs[n] = -weights[n];
			n++;
		}
	}
	MinSpanTree(tree, nClusters, m, edges, costs);

	nSeperators = nClusters - 1;
	seperatorEdges = (int **) R_allocArray<int>(nSeperators, 2);
	seperatorSize = (int *) R_allocVector<int>(nSeperators);
	n = 0;
	for (int i = 0; i < m; i++)
	{
		if (tree[i])
		{
			seperatorEdges[n][0] = edges[i];
			seperatorEdges[n][1] = edges[i + m];
			seperatorSize[n] = weights[i];
			n++;
		}
	}
	C_freeVector(overlap);
	C_freeVector(tree);
	C_freeVector(edges);
	C_freeVector(weights);
	C_freeVector(costs);

	seperators = (int **) R_allocArray2<int>(nSeperators, seperatorSize);
	for (int i = 0; i < nSeperators; i++)
	{
		n1 = seperatorEdges[i][0];
		n2 = seperatorEdges[i][1];
		Intersection(seperators[i], clusters[n1], clusterSize[n1], clusters[n2], clusterSize[n2]);
	}

	nClusterStates = (int *) R_allocVector<int>(nClusters);
	for (int i = 0; i < nClusters; i++)
	{
		nClusterStates[i] = 1;
		for (int j = 0; j < clusterSize[i]; j++)
			nClusterStates[i] *= original.nStates[clusters[i][j]];
	}

	nSeperatorStates = (int *) R_allocVector<int>(nSeperators);
	for (int i = 0; i < nSeperators; i++)
	{
		nSeperatorStates[i] = 1;
		for (int j = 0; j < seperatorSize[i]; j++)
			nSeperatorStates[i] *= original.nStates[seperators[i][j]];
	}

	clusterBel = (double **) R_allocArray2<double>(nClusters, nClusterStates);
	seperatorBel = (double **) R_allocArray2<double>(nSeperators, nSeperatorStates);
}
