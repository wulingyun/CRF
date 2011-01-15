#include "CRF.h"

SEXP Decode_Tree(SEXP _crf)
{
	CRFinfo crf;
	openCRF(&crf, _crf);

	SEXP _labels;
	PROTECT(_labels = NEW_INTEGER(crf.nNodes));
	int *labels = INTEGER_POINTER(_labels);
	setValues(_labels, labels, -1);

	_Decode_Tree(&crf, labels);

	UNPROTECT(1);
	closeCRF(&crf);

	return(_labels);
}

void _Decode_Tree(CRFinfo *crf, int *labels)
{
	/* Tree BP */

	double *messages_1 = (double *) R_alloc(crf->maxState * crf->nEdges, sizeof(double)); // Messages from n2 to n1 at edge (n1, n2)
	double *messages_2 = (double *) R_alloc(crf->maxState * crf->nEdges, sizeof(double)); // Messages from n1 to n2 at edge (n1, n2)
	TreeBP_max(crf, messages_1, messages_2);

	/* Node beliefs */

	double *nodeBel = (double *) R_alloc(crf->nNodes * crf->maxState, sizeof(double));
	Message2NodeBelief(crf, messages_1, messages_2, nodeBel);

	/* Labels */

	double maxBel, *p_nodeBel;
	for (int i = 0; i < crf->nNodes; i++)
	{
		maxBel = -1;
		p_nodeBel = nodeBel + i;
		for (int j = 0; j < crf->nStates[i]; j++)
		{
			if (p_nodeBel[0] > maxBel)
			{
				maxBel = p_nodeBel[0];
				labels[i] = j;
			}
			p_nodeBel += crf->nNodes;
		}
	}

	for (int i = 0; i < crf->nNodes; i++)
		labels[i]++;
}
