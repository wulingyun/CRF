#include "CRF.h"

/* Junction tree init */

void CRF::JunctionTreeInit()
{
	int dim[] = {nNodes, nNodes};
	int **adj = (int **) allocArray<int, 2>(dim);
	int **neighbors = (int **) allocArray<int, 2>(dim);
	int **clusters = (int **) allocArray<int, 2>(dim);
	int *nNeighbors = (int *) R_alloc(nNodes, sizeof(int));
	int *clusterSize = (int *) R_alloc(nNodes, sizeof(int));
	int nClusters = 0;
	for (int i = 0; i < nNodes; i++)
	{
		for (int j = 0; j < nNodes; j++)
		{
			adj[i][j] = 0;
			neighbors[i][j] = 0;
			clusters[i][j] = 0;
		}
		nNeighbors[i] = 0;
		clusterSize[i] = 0;
	}

	int n1, n2;
	for (int i = 0; i < nEdges; i++)
	{
		n1 = EdgesBegin(i);
		n2 = EdgesEnd(i);
		adj[n1][n2] = adj[n2][n1] = 1;
		neighbors[n1][nNeighbors[n1]++] = n2;
		neighbors[n2][nNeighbors[n2]++] = n1;
	}

	int n, m;
	while (1)
	{
		n = -1;
		m = nNodes;
		for (int i = 0; i < nNodes; i++)
		{
			if (nNeighbors[i] > 0 && nNeighbors[i] < m)
			{
				n = i;
				m = nNeighbors[i];
			}
		}
		if (n < 0)
			break;

		for (int i = 0; i < m; i++)
		{
			n1 = neighbors[n][i];
			for (int j = i+1; j < m; j++)
			{
				n2 = neighbors[n][j];
				if (adj[n1][n2] == 0)
				{
					adj[n1][n2] = adj[n2][n1] = 1;
					neighbors[n1][nNeighbors[n1]++] = n2;
					neighbors[n2][nNeighbors[n2]++] = n1;
				}
			}
			for (int j = 0; j < nNeighbors[n1]; j++)
			{
				if (neighbors[n1][j] == n)
				{
					for (int k = j; k < nNeighbors[n1]-1; k++)
					{
						neighbors[n1][k] = neighbors[n1][k+1];
					}
					nNeighbors[n1]--;
					break;
				}
			}
			adj[n1][n] = adj[n][n1] = 0;
			clusters[nClusters][i] = neighbors[n][i];
		}
		nNeighbors[n] = 0;
		clusters[nClusters][m] = n;
		clusterSize[nClusters++] = m+1;
	}
	int treeWidth = 0;
	for (int i = 0; i < nClusters; i++)
	{
		if (clusterSize[i] > treeWidth)
			treeWidth = clusterSize[i];
	}
	treeWidth--;
	Rprintf("%d, %d\n", nClusters, treeWidth);
}
