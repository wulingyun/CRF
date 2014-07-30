#include <R.h>
#include <Rdefines.h>

/* macros */

#define min(a, b) a < b ? a : b;
#define max(a, b) a > b ? a : b;

/* initialize the list */

template <class T>
inline void SetValues(SEXP r, T *c, T v)
{
	for (int i = 0; i < length(r); i++)
		c[i] = v;
};

/* get the variable */

inline SEXP GetVar(SEXP env, const char *name)
{
	return findVar(install(name), env);
};

/* get the variable */

inline void SetVar(SEXP env, const char *name, SEXP var)
{
	defineVar(install(name), var, env);
};

/* get/set the list element */

inline SEXP GetListElement(SEXP list, int i)
{
  if (i >= 0 && i < length(list)) return VECTOR_ELT(list, i);
  return R_NilValue;
}

SEXP GetListElement(SEXP list, const char *tag);
void SetListElement(SEXP list, int i, const char *tag, SEXP value);

/* set dim of array */

void SetDim2(SEXP array, int x1, int x2);
void SetDim3(SEXP array, int x1, int x2, int x3);

/* sample from discrete distribution */

int SampleFrom(int n, double *prob);

/* minimum weight spanning tree using Kruskal algorithm */

int MinSpanTree(int *tree, int nNodes, int nEdges, int *edges, double *costs, int node_index_from = 1);

/* utils for ascending ordered vector */

int Intersection(int *overlap, int *vector1, int size1, int *vector2, int size2);
void Insert(int *vector, int &size, int v);
void Remove(int *vector, int &size, int v);

/* swap variables */

template <class T>
inline void swap(T &a, T &b)
{
	T temp = a;
	a = b;
	b = temp;
};

/* allocate vector */

template <class T>
inline T *R_allocVector(int n)
{
	T *vector = (T *) R_alloc(n, sizeof(T));
	return vector;
};

/* allocate vector (Calloc version) */

template <class T>
inline T *C_allocVector(int n)
{
	T *vector = (T *) Calloc(n, T);
	return vector;
};

template <class T>
inline void C_freeVector(T *vector)
{
	Free(vector);
};

/* allocate multidimensional array */

template <class T, int n>
inline T **R_allocArray(int dim[n])
{
	int size1, size2;
	T **array, **sub1, **sub2, **tmp;
	size1 = dim[0];
	array = sub1 = (T **) R_alloc(size1, sizeof(T *));
	for (int i = 1; i < n-1; i++)
	{
		size2 = size1 * dim[i];
		tmp = sub2 = (T **) R_alloc(size2, sizeof(T *));
		for (int j = 0; j < size1; j++)
		{
			sub1[j] = (T *) tmp;
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
inline T **R_allocArray(int dim1, int dim2)
{
	int dim[] = {dim1, dim2};
	return R_allocArray<T, 2>(dim);
};

template <class T>
inline T **R_allocArray(int dim1, int dim2, int dim3)
{
	int dim[] = {dim1, dim2, dim3};
	return R_allocArray<T, 3>(dim);
};

/* allocate multidimensional array (Calloc version) */

template <class T, int n>
inline T **C_allocArray(int dim[n])
{
	int size1, size2;
	T **array, **sub1, **sub2, **tmp;
	size1 = dim[0];
	array = sub1 = (T **) Calloc(size1, T *);
	for (int i = 1; i < n-1; i++)
	{
		size2 = size1 * dim[i];
		tmp = sub2 = (T **) Calloc(size2, T *);
		for (int j = 0; j < size1; j++)
		{
			sub1[j] = (T *) tmp;
			tmp += dim[i];
		}
		size1 = size2;
		sub1 = sub2;
	}
	T *block = (T *) Calloc(size1 * dim[n-1], T);
	for (int i = 0; i < size1; i++)
	{
		sub1[i] = block;
		block += dim[n-1];
	}
	return array;
};

template <class T>
inline T **C_allocArray(int dim1, int dim2)
{
	int dim[] = {dim1, dim2};
	return C_allocArray<T, 2>(dim);
};

template <class T>
inline T **C_allocArray(int dim1, int dim2, int dim3)
{
	int dim[] = {dim1, dim2, dim3};
	return C_allocArray<T, 3>(dim);
};

template <class T, int n>
inline void C_freeArray(T **array)
{
	void *p;
	for (int i = 0; i < n-1; i++)
	{
		p = array;
		array = (T **) array[0];
		Free(p);
	}
  Free(array);
};

/* allocate 2D array with varied dim2 */

template <class T>
inline T **R_allocArray2(int dim1, int *dim2)
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

/* allocate 2D array with varied dim2 (Calloc version) */

template <class T>
inline T **C_allocArray2(int dim1, int *dim2)
{
	T *block, **array;
	int array_size;
	array_size = 0;
	for (int i = 0; i < dim1; i++)
		array_size += dim2[i];
	block = (T *) Calloc(array_size, T);
	array = (T **) Calloc(dim1, T *);
	for (int i = 0; i < dim1; i++)
	{
		array[i] = block;
		block += dim2[i];
	}
	return array;
};

template <class T>
inline void C_freeArray2(T **array)
{
	Free(array[0]);
	Free(array);
};
