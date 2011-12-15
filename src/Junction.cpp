#include "CRF.h"

/* Get intersection of two ascending ordered vectors */

int Intersection(int *overlap, int *vector1, int size1, int *vector2, int size2)
{
	int n, i1, i2;
	n = i1 = i2 = 0;
	while (i1 < size1 && i2 < size2)
	{
		if (vector1[i1] == vector2[i2])
		{
			overlap[n++] = vector1[i1++];
			i2++;
		}
		else if (vector1[i1] < vector2[i2])
			i1++;
		else
			i2++;
	}
	return n;
}

/* Insert element to ascending ordered vector */

void Insert(int *vector, int &size, int v)
{
	int k = size;
	for (int i = 0; i < size; i++)
	{
		if (vector[i] > v)
		{
			for (int j = size; j > i; j--)
				vector[j] = vector[j-1];
			k = i;
			break;
		}
	}
	vector[k] = v;
	size++;
}

/* Remove element from ascending ordered vector */

void Remove(int *vector, int &size, int v)
{

	for (int i = 0; i < size; i++)
	{
		if (vector[i] == v)
		{
			for (int j = i; j < size-1; j++)
				vector[j] = vector[j+1];
			size--;
			break;
		}
	}
}

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
		Insert(neighbors[n1], nNeighbors[n1], n2);
		Insert(neighbors[n2], nNeighbors[n2], n1);
		adj[n1][n2] = adj[n2][n1] = 1;
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

	int n, m, maxMissingEdges = nNodes * nNodes;
	int treeWidth = 0;
	int *overlap = (int *) R_alloc(nNodes, sizeof(int));
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
			Insert(clusters[nClusters], clusterSize[nClusters], n1);
		}
		nMissingEdges[n] = -1;
		nNeighbors[n] = 0;
		Insert(clusters[nClusters], clusterSize[nClusters], n);
		if (clusterSize[nClusters] > treeWidth)
			treeWidth = clusterSize[nClusters];
		nClusters++;
	}
	treeWidth--;
	Rprintf("%d, %d\n", nClusters, treeWidth);
}
