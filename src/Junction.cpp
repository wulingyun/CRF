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

	int *nMissingEdges = (int *) R_alloc(nNodes, sizeof(int));
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

	int n, m;
	int treeWidth = 0;
	while (1)
	{
		n = -1;
		m = nNodes * nNodes;
		for (int i = 0; i < nNodes; i++)
		{
			if (nNeighbors[i] > 0 && nMissingEdges[i] < m)
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
					m = 0;
					for (int k1 = 0; k1 < nNeighbors[n1]; k1++)
					{
						for (int k2 = 0; k2 < nNeighbors[n2]; k2++)
						{
							if (neighbors[n1][k1] == neighbors[n2][k2])
							{
								m++;
								nMissingEdges[neighbors[n1][k1]]--;
								break;
							}
						}
					}
					nMissingEdges[n1] += nNeighbors[n1]-m;
					nMissingEdges[n2] += nNeighbors[n2]-m;
					neighbors[n1][nNeighbors[n1]++] = n2;
					neighbors[n2][nNeighbors[n2]++] = n1;
					adj[n1][n2] = adj[n2][n1] = 1;
				}
			}
		}
		for (int j1 = 0; j1 < nNeighbors[n]; j1++)
		{
			n1 = neighbors[n][j1];
			clusters[nClusters][j1] = n1;
			m = 0;
			for (int k1 = 0; k1 < nNeighbors[n1]; k1++)
			{
				for (int k = 0; k < nNeighbors[n]; k++)
				{
					if (neighbors[n1][k1] == neighbors[n][k])
					{
						m++;
						break;
					}
				}
			}
			nMissingEdges[n1] -= nNeighbors[n1]-1-m;
			for (int k1 = 0; k1 < nNeighbors[n1]; k1++)
			{
				if (neighbors[n1][k1] == n)
				{
					for (int j = k1; j < nNeighbors[n1]-1; j++)
						neighbors[n1][j] = neighbors[n1][j+1];
					nNeighbors[n1]--;
					break;
				}
			}
		}
		clusters[nClusters][nNeighbors[n]] = n;
		clusterSize[nClusters++] = nNeighbors[n]+1;
		nNeighbors[n] = 0;

		if (clusterSize[nClusters-1] > treeWidth)
			treeWidth = clusterSize[nClusters-1];
		//Rprintf("%d, %d\n", nClusters, treeWidth);
	}
	treeWidth--;
	Rprintf("%d\n", treeWidth);
}
