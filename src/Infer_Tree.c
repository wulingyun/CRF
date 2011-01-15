#include "CRF.h"

SEXP Infer_Tree(SEXP _crf)
{
	CRFinfo crf;
	openCRF(&crf, _crf);

	SEXP _nodeBel, _edgeBel, _logZ;
	PROTECT(_nodeBel = NEW_NUMERIC(crf.nNodes * crf.maxState));
	PROTECT(_edgeBel = NEW_NUMERIC(crf.maxState * crf.maxState * crf.nEdges));
	PROTECT(_logZ = NEW_NUMERIC(1));
	setDim2(_nodeBel, crf.nNodes, crf.maxState);
	setDim3(_edgeBel, crf.maxState, crf.maxState, crf.nEdges);
	double *nodeBel = NUMERIC_POINTER(_nodeBel);
	double *edgeBel = NUMERIC_POINTER(_edgeBel);
	double *logZ = NUMERIC_POINTER(_logZ);
	setValues(_nodeBel, nodeBel, 0);
	setValues(_edgeBel, edgeBel, 0);
	*logZ = 0;

	_Infer_Tree(&crf, nodeBel, edgeBel, logZ);

	SEXP _belief;
	PROTECT(_belief = NEW_LIST(3));
	setListElement(_belief, 0, "node.bel", _nodeBel);
	setListElement(_belief, 1, "edge.bel", _edgeBel);
	setListElement(_belief, 2, "logZ", _logZ);

	UNPROTECT(4);
	closeCRF(&crf);

	return(_belief);
}

void _Infer_Tree(CRFinfo *crf, double *nodeBel, double *edgeBel, double *logZ)
{
	/* Tree BP */

	double *messages_1 = (double *) R_alloc(crf->maxState * crf->nEdges, sizeof(double)); // Messages from n2 to n1 at edge (n1, n2)
	double *messages_2 = (double *) R_alloc(crf->maxState * crf->nEdges, sizeof(double)); // Messages from n1 to n2 at edge (n1, n2)
	TreeBP(crf, messages_1, messages_2);

	/* Beliefs */

	Message2NodeBelief(crf, messages_1, messages_2, nodeBel);
	Message2EdgeBelief(crf, messages_1, messages_2, nodeBel, edgeBel);

	/* Bethe free energy */

	double nodeEnergy, nodeEntropy, edgeEnergy, edgeEntropy;
	nodeEnergy = nodeEntropy = edgeEnergy = edgeEntropy = 0;

	double entropy;
	double *p_nodeBel, *p_nodePot;
	for (int i = 0; i < crf->nNodes; i++)
	{
		entropy = 0;
		p_nodeBel = nodeBel + i;
		p_nodePot = crf->nodePot + i;
		for (int j = 0; j < crf->nStates[i]; j++)
		{
			if (p_nodeBel[0] > 0)
			{
				nodeEnergy -= p_nodeBel[0] * log(p_nodePot[0]);
				entropy += p_nodeBel[0] * log(p_nodeBel[0]);
			}
			p_nodeBel += crf->nNodes;
			p_nodePot += crf->nNodes;
		}
		nodeEntropy += (crf->nAdj[i] - 1) * entropy;
	}

	int n1, n2;
	double *p_edgeBel, *p0_edgeBel = edgeBel;
	double *p_edgePot, *p0_edgePot = crf->edgePot;
	for (int i = 0; i < crf->nEdges; i++)
	{
		n1 = crf->edges[i] - 1;
		n2 = crf->edges[i + crf->nEdges] - 1;
		p_edgeBel = p0_edgeBel;
		p_edgePot = p0_edgePot;
		for (int j = 0; j < crf->nStates[n2]; j++)
		{
			for (int k = 0; k < crf->nStates[n1]; k++)
			{
				if (p_edgeBel[k] > 0)
				{
					edgeEnergy -= p_edgeBel[k] * log(p_edgePot[k]);
					edgeEntropy -= p_edgeBel[k] * log(p_edgeBel[k]);
				}
			}
			p_edgeBel += crf->maxState;
			p_edgePot += crf->maxState;
		}
		p0_edgeBel += crf->maxState * crf->maxState;
		p0_edgePot += crf->maxState * crf->maxState;
	}

	*logZ = - nodeEnergy + nodeEntropy - edgeEnergy + edgeEntropy;
}
