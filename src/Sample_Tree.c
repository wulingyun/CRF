#include "CRF.h"

SEXP Sample_Tree(SEXP _crf, SEXP _size)
{
	CRFinfo crf;
	openCRF(&crf, _crf);

	PROTECT(_size = AS_INTEGER(_size));
	int size = INTEGER_POINTER(_size)[0];

	SEXP _samples;
	PROTECT(_samples = NEW_INTEGER(size * crf.nNodes));
	setDim2(_samples, size, crf.nNodes);
	int *samples = INTEGER_POINTER(_samples);
	setValues(_samples, samples, 0);

	_Sample_Tree(&crf, size, samples);

	UNPROTECT(2);
	closeCRF(&crf);

	return(_samples);
}

void _Sample_Tree(CRFinfo *crf, int size, int *samples)
{
	int *y = (int *) R_alloc(crf->nNodes, sizeof(int));
	for (int i = 0; i < crf->nNodes; i++)
		y[i] = 0;

	/* Tree BP */

	double *messages_1 = (double *) R_alloc(crf->maxState * crf->nEdges, sizeof(double)); // Messages from n2 to n1 at edge (n1, n2)
	double *messages_2 = (double *) R_alloc(crf->maxState * crf->nEdges, sizeof(double)); // Messages from n1 to n2 at edge (n1, n2)
	TreeBP(crf, messages_1, messages_2);

	/* Beliefs */

	double *nodeBel = (double *) R_alloc(crf->nNodes * crf->maxState, sizeof(double));
	double *edgeBel = (double *) R_alloc(crf->maxState * crf->maxState * crf->nEdges, sizeof(double));
	Message2NodeBelief(crf, messages_1, messages_2, nodeBel);
	Message2EdgeBelief(crf, messages_1, messages_2, nodeBel, edgeBel);

	/* Sampling */

	int nOrdered = 0;
	int *ordered = (int *) R_alloc(crf->nNodes, sizeof(int));
	int *order = (int *) R_alloc(crf->nNodes, sizeof(int));
	int *parentEdge = (int *) R_alloc(crf->nNodes, sizeof(int));

	int nQueue;
	int *queue = (int *) R_alloc(crf->nNodes, sizeof(int *));
	int s, e, n;

	for (int i = 0; i < crf->nNodes; i++)
		ordered[i] = queue[i] = 0;

	int n1, n2;
	for (int i = 0; i < crf->nNodes; i++)
	{
		if (ordered[i])
			continue;

		ordered[i] = 1;
		order[nOrdered] = i;
		parentEdge[nOrdered] = -1;
		nOrdered++;

		queue[i] = 1;
		nQueue = 1;
		while (nQueue > 0)
		{
			n1 = 0;
			for (int j = 0; j < crf->nNodes; j++)
				if (queue[j])
				{
					n1 = j;
					queue[j] = 0;
					nQueue--;
					break;
				}

			for (int j = 0; j < crf->nAdj[n1]; j++)
			{
				n2 = crf->adjNodes[n1][j] - 1;
				if (ordered[n2])
					continue;

				ordered[n2] = 1;
				order[nOrdered] = n2;
				parentEdge[nOrdered] = crf->adjEdges[n1][j] - 1;
				nOrdered++;

				queue[n2] = 1;
				nQueue++;
			}
		}
	}

	double sumProb, *prob = (double *) R_alloc(crf->maxState, sizeof(double));
	double *p_nodeBel, *p_edgeBel;

	GetRNGstate();
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < crf->nNodes; j++)
		{
			n = order[j];
			e = parentEdge[j];
			if (e == -1)
			{
				p_nodeBel = nodeBel + n;
				for (int k = 0; k < crf->nStates[n]; k++)
				{
					prob[k] = p_nodeBel[0];
					p_nodeBel += crf->nNodes;
				}
			}
			else
			{
				sumProb = 0;
				if (crf->edges[e] - 1 == n)
				{
					s = crf->edges[e + crf->nEdges] - 1;
					p_edgeBel = edgeBel + crf->maxState * (y[s] + crf->maxState * e);
					for (int k = 0; k < crf->nStates[n]; k++)
					{
						prob[k] = p_edgeBel[k];
						sumProb += prob[k];
					}
				}
				else
				{
					s = crf->edges[e] - 1;
					p_edgeBel = edgeBel + y[s] + crf->maxState * crf->maxState * e;
					for (int k = 0; k < crf->nStates[n]; k++)
					{
						prob[k] = p_edgeBel[0];
						sumProb += prob[k];
						p_edgeBel += crf->maxState;
					}
				}
				for (int k = 0; k < crf->nStates[n]; k++)
					prob[k] /= sumProb;
			}
			y[n] = sample(crf->nStates[n], prob);
		}

		for (int j = 0; j < crf->nNodes; j++)
			samples[i + size * j] = y[j] + 1;
	}
	PutRNGstate();
}
