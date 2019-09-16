#include <R.h>
#include <Rdefines.h>

/* sample from discrete distribution */

int SampleFrom(int n, double *prob);

/* minimum weight spanning tree using Kruskal algorithm */

int MinSpanTree(int *tree, int nNodes, int nEdges, int *edges, double *costs, int node_index_from = 1);

/* utils for ascending ordered vector */

int Intersection(int *overlap, int *vector1, int size1, int *vector2, int size2);
void Insert(int *vector, int &size, int v);
void Remove(int *vector, int &size, int v);
