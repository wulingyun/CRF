#include "CRF.h"

/* sample from discret distribution */

int SampleFrom(int n, double *prob)
{
	int select = n-1;
	double cutoff = unif_rand();
	double cumulativeProb = 0;
	for (int i = 0; i < n; i++)
	{
		cumulativeProb += prob[i];
		if (cumulativeProb > cutoff)
		{
			select = i;
			break;
		}
	}
	return select;
}

/* minimum weight spanning tree using Kruskal algorithm */

int MinSpanTree(int *tree, int nNodes, int nEdges, int *edges, double *costs, int node_index_from)
{
	int *index = (int *) C_allocVector<int>(nEdges);
	for (int i = 0; i < nEdges; i++)
	{
		tree[i] = 0;
		index[i] = i;
	}
	rsort_with_index(costs, index, nEdges);

	int *label = (int *) C_allocVector<int>(nNodes);
	for (int i = 0; i < nNodes; i++)
		label[i] = i;

	int n = 0, n1, n2;
	for (int i = 0; i < nEdges; i++)
	{
		n1 = label[edges[index[i]] - node_index_from];
		n2 = label[edges[index[i] + nEdges] - node_index_from];
		if (n1 != n2)
		{
			for (int j = 0; j < nNodes; j++)
				if (label[j] == n2)
					label[j] = n1;
			tree[index[i]] = 1;
			if (++n >= nNodes - 1)
				break;
		}
	}

	C_freeVector(index);
	C_freeVector(label);

	return n;
}

/* Get intersection of two ascending ordered vectors */

int Intersection(int *overlap, int *vector1, int size1, int *vector2, int size2)
{
	int n, i1, i2;
	if (vector1[0] > vector2[size2-1] || vector2[0] > vector1[size1-1])
		return 0;
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

/* Calculate Frequency */

SEXP Calc_Frequency(SEXP _v, SEXP _n)
{
	PROTECT(_v = AS_INTEGER(_v));
	PROTECT(_n = AS_INTEGER(_n));
	int *v = INTEGER_POINTER(_v);
	int n = INTEGER_POINTER(_n)[0];
	int m = length(_v);

	SEXP _freq;
	PROTECT(_freq = NEW_INTEGER(n));
	int *freq = INTEGER_POINTER(_freq);
	SetValues(_freq, freq, 0);

	for (int i = 0; i < m; i++)
		if (v[i] > 0 && v[i] <= n)
			freq[v[i]-1]++;

	UNPROTECT(3);
	return(_freq);
}
