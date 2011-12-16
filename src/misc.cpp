#include "CRF.h"

/* get the list element named tag, or return NULL */

SEXP GetListElement(SEXP list, const char *tag)
{
	SEXP value = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	for (int i = 0; i < length(list); i++)
		if (strcmp(CHAR(STRING_ELT(names, i)), tag) == 0)
		{
			value = VECTOR_ELT(list, i);
			break;
		}
	return value;
}

/* set the list element */

void SetListElement(SEXP list, int i, const char *tag, SEXP value)
{
	SEXP _names = getAttrib(list, R_NamesSymbol);
	if (_names == R_NilValue)
	{
		PROTECT(_names = NEW_STRING(length(list)));
		SET_STRING_ELT(_names, i, mkChar(tag));
		setAttrib(list, R_NamesSymbol, _names);
		UNPROTECT_PTR(_names);
	}
	else
		SET_STRING_ELT(_names, i, mkChar(tag));
	SET_VECTOR_ELT(list, i, value);
}

/* set dim of array */

void SetDim2(SEXP array, int x1, int x2)
{
	SEXP _dim;
	PROTECT(_dim = NEW_INTEGER(2));
	INTEGER_POINTER(_dim)[0] = x1;
	INTEGER_POINTER(_dim)[1] = x2;
	SET_DIM(array, _dim);
	UNPROTECT_PTR(_dim);
}

void SetDim3(SEXP array, int x1, int x2, int x3)
{
	SEXP _dim;
	PROTECT(_dim = NEW_INTEGER(3));
	INTEGER_POINTER(_dim)[0] = x1;
	INTEGER_POINTER(_dim)[1] = x2;
	INTEGER_POINTER(_dim)[2] = x3;
	SET_DIM(array, _dim);
	UNPROTECT_PTR(_dim);
}

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
			for (int j = 0; j < nNodes; j++)
				if (label[j] == label[n2])
					label[j] = label[n1];
			tree[index[i]] = 1;
			if (++n >= nNodes - 1)
				break;
		}
	}
}
