#include "CRF.h"

SEXP Sample_Tree(SEXP _crf, SEXP _size)
{
	CRF crf(_crf);
	crf.Init_Samples(_size);
	crf.Init_NodeBel();
	crf.Init_EdgeBel();
	crf.Sample_Tree();
	return(crf._samples);
}

void CRF::Sample_Tree(int size)
{
	void *vmax = vmaxget(); 

	if (size <= 0)
		size = nSamples;
	else if (size > nSamples)
		Init_Samples(size);

	int *y = (int *) R_alloc(nNodes, sizeof(int));
	for (int i = 0; i < nNodes; i++)
		y[i] = 0;

	/* Tree BP */

	TreeBP();
	Messages2EdgeBel();

	/* Sampling */

	int nOrdered = 0, nQueue = 0;
	int *ordered = (int *) R_alloc(nNodes, sizeof(int));
	int *order = (int *) R_alloc(nNodes, sizeof(int));
	int *parentEdge = (int *) R_alloc(nNodes, sizeof(int));
	int *queue = (int *) R_alloc(nNodes, sizeof(int));

	for (int i = 0; i < nNodes; i++)
		ordered[i] = 0;

	int s, e, n, n1, n2;
	for (int i = 0; i < nNodes; i++)
	{
		if (ordered[i])
			continue;

		ordered[i] = 1;
		order[nOrdered] = i;
		parentEdge[nOrdered] = -1;
		nOrdered++;

		queue[nQueue++] = i;
		while (nQueue > 0)
		{
			n1 = queue[--nQueue];
			for (int j = 0; j < nAdj[n1]; j++)
			{
				n2 = AdjNodes(n1, j);
				if (ordered[n2])
					continue;

				ordered[n2] = 1;
				order[nOrdered] = n2;
				parentEdge[nOrdered] = AdjEdges(n1, j);
				nOrdered++;

				queue[nQueue++] = n2;
			}
		}
	}

	double sumProb, *prob = (double *) R_alloc(maxState, sizeof(double));

	GetRNGstate();
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < nNodes; j++)
		{
			n = order[j];
			e = parentEdge[j];
			if (e == -1)
			{
				for (int k = 0; k < nStates[n]; k++)
					prob[k] = NodeBel(n, k);
			}
			else
			{
				sumProb = 0;
				if (EdgesBegin(e) == n)
				{
					s = EdgesEnd(e);
					for (int k = 0; k < nStates[n]; k++)
					{
						prob[k] = EdgeBel(e, k, y[s]);
						sumProb += prob[k];
					}
				}
				else
				{
					s = EdgesBegin(e);
					for (int k = 0; k < nStates[n]; k++)
					{
						prob[k] = EdgeBel(e, y[s], k);
						sumProb += prob[k];
					}
				}
				for (int k = 0; k < nStates[n]; k++)
					prob[k] /= sumProb;
			}
			y[n] = SampleFrom(nStates[n], prob);
		}

		for (int j = 0; j < nNodes; j++)
			Samples(i, j) = y[j] + 1;
	}
	PutRNGstate();

	vmaxset(vmax); 
}
