#include "CRF.h"

SEXP Decode_Tree(SEXP _crf)
{
	CRF crf(_crf);

	SEXP _labels;
	PROTECT(_labels = NEW_INTEGER(crf.nNodes));
	int *labels = INTEGER_POINTER(_labels);
	setValues(_labels, labels, -1);

	crf.Decode_Tree(labels);

	UNPROTECT(1);

	return(_labels);
}

void CRF::Decode_Tree(int *labels)
{
	/* Tree BP */

	double *messages_1 = (double *) R_alloc(maxState * nEdges, sizeof(double)); // Messages from n2 to n1 at edge (n1, n2)
	double *messages_2 = (double *) R_alloc(maxState * nEdges, sizeof(double)); // Messages from n1 to n2 at edge (n1, n2)
	TreeBP_max(messages_1, messages_2);

	/* Node beliefs */

	double *nodeBel = (double *) R_alloc(nNodes * maxState, sizeof(double));
	Message2NodeBelief(messages_1, messages_2, nodeBel);

	/* Labels */

	double maxBel, *p_nodeBel;
	for (int i = 0; i < nNodes; i++)
	{
		maxBel = -1;
		p_nodeBel = nodeBel + i;
		for (int j = 0; j < nStates[i]; j++)
		{
			if (p_nodeBel[0] > maxBel)
			{
				maxBel = p_nodeBel[0];
				labels[i] = j;
			}
			p_nodeBel += nNodes;
		}
	}

	for (int i = 0; i < nNodes; i++)
		labels[i]++;
}
