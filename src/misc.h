#include <R.h>
#include <Rdefines.h>

/* initialize the list */
template <class T>
void SetValues(SEXP r, T *c, T v)
{
	for (int i = 0; i < length(r); i++)
		c[i] = v;
};

/* get/set the list element */
SEXP GetListElement(SEXP list, const char *tag);
void SetListElement(SEXP list, int i, const char *tag, SEXP value);

/* set dim of array */
void SetDim2(SEXP array, int x1, int x2);
void SetDim3(SEXP array, int x1, int x2, int x3);

/* sample from discrete distribution */
int SampleFrom(int n, double *prob);

/* minimum weight spanning tree using Kruskal algorithm */
void MinSpanTree(int *tree, int nNodes, int nEdges, int *edges, double *costs);

/* swap variables */
template <class T>
void swap(T &a, T &b)
{
	T temp = a;
	a = b;
	b = temp;
};

/* allocate multidimensional array */
template <class T, int n>
void **allocArray(int dim[n])
{
	int size1, size2;
	void **array, **sub1, **sub2, **tmp;
	size1 = dim[0];
	array = sub1 = (void **) R_alloc(size1, sizeof(void *));
	for (int i = 1; i < n-1; i++)
	{
		size2 = size1 * dim[i];
		tmp = sub2 = (void **) R_alloc(size2, sizeof(void *));
		for (int j = 0; j < size1; j++)
		{
			sub1[j] = (void *) tmp;
			tmp += dim[i];
		}
		size1 = size2;
		sub1 = sub2;
	}
	T *block = (T *) R_alloc(size1 * dim[n-1], sizeof(T));
	for (int i = 0; i < size1; i++)
	{
		sub1[i] = block;
		block += dim[n-1];
	}
	return array;
};

template <class T>
T **allocArray2(int dim1, int *dim2)
{
	T *block, **array;
	int array_size;
	array_size = 0;
	for (int i = 0; i < dim1; i++)
		array_size += dim2[i];
	block = (T *) R_alloc(array_size, sizeof(T));
	array = (T **) R_alloc(dim1, sizeof(T *));
	for (int i = 0; i < dim1; i++)
	{
		array[i] = block;
		block += dim2[i];
	}
	return array;
};
