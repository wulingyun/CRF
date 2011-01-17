#include "CRF.h"

SEXP Infer_Tree(SEXP _crf)
{
	CRF crf(_crf);

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

	crf.Infer_Tree(nodeBel, edgeBel, logZ);

	SEXP _belief;
	PROTECT(_belief = NEW_LIST(3));
	setListElement(_belief, 0, "node.bel", _nodeBel);
	setListElement(_belief, 1, "edge.bel", _edgeBel);
	setListElement(_belief, 2, "logZ", _logZ);

	UNPROTECT(4);

	return(_belief);
}

void CRF::Infer_Tree(double *nodeBel, double *edgeBel, double *logZ)
{
	double *messages_1 = (double *) R_alloc(maxState * nEdges, sizeof(double)); // Messages from n2 to n1 at edge (n1, n2)
	double *messages_2 = (double *) R_alloc(maxState * nEdges, sizeof(double)); // Messages from n1 to n2 at edge (n1, n2)
	TreeBP(messages_1, messages_2);
	Message2NodeBelief(messages_1, messages_2, nodeBel);
	Message2EdgeBelief(messages_1, messages_2, nodeBel, edgeBel);
	BetheFreeEnergy(nodeBel, edgeBel, logZ);
}
