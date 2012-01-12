#include "CRF.h"

SEXP Update_Pot(SEXP _crf, SEXP _nf, SEXP _ef)
{
	CRF crf(_crf);
	double *nf = NUMERIC_POINTER(AS_NUMERIC(_nf));
	double *ef = NUMERIC_POINTER(AS_NUMERIC(_ef));

	int nNodeFea = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.nf")))[0];
	int nEdgeFea = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.ef")))[0];

	SEXP _par;
	double *par;
	PROTECT(_par = AS_NUMERIC(GetListElement(_crf, "par")));
	par = NUMERIC_POINTER(_par);

	SEXP _nodePar, _edgePar, _temp;
	int *nodePar, **edgePar;
	PROTECT(_nodePar = AS_INTEGER(GetListElement(_crf, "node.par")));
	PROTECT(_edgePar = AS_LIST(GetListElement(_crf, "edge.par")));
	nodePar = INTEGER_POINTER(_nodePar);
	edgePar = (int **) R_alloc(crf.nEdges, sizeof(int *));
	for (int i = 0; i < crf.nEdges; i++)
	{
		PROTECT(_temp = AS_INTEGER(VECTOR_ELT(_edgePar, i)));
		edgePar[i] = INTEGER_POINTER(_temp);
	}

	for (int i = 0; i < crf.nNodes * crf.maxState; i++)
		crf.nodePot[i] = 0;
	for (int i = 0; i < crf.nEdges; i++)
		for (int j = 0; j < crf.nEdgeStates[i]; j++)
			crf.edgePot[i][j] = 0;

	for (int i = 0; i < crf.nNodes; i++)
		for (int j = 0; j < nNodeFea; j++)
		{
			double f = nf[j + nNodeFea * i];
			if (f != 0)
				for (int k = 0; k < crf.nStates[i]; k++)
				{
					int p = nodePar[i + crf.nNodes * (k + crf.maxState * j)] - 1;
					if (p >= 0)
						crf.nodePot[i + crf.nNodes * k] += f * par[p];
				}
		}

	for (int i = 0; i < crf.nEdges; i++)
		for (int j = 0; j < nEdgeFea; j++)
		{
			double f = ef[j + nEdgeFea * i];
			if (f != 0)
				for (int k = 0; k < crf.nEdgeStates[i]; k++)
				{
					int p = edgePar[i][k + crf.nEdgeStates[i] * j] - 1;
					if (p >= 0)
						crf.edgePot[i][k] += f * par[p];
				}
		}

	for (int i = 0; i < crf.nNodes * crf.maxState; i++)
		crf.nodePot[i] = exp(crf.nodePot[i]);
	for (int i = 0; i < crf.nEdges; i++)
		for (int j = 0; j < crf.nEdgeStates[i]; j++)
			crf.edgePot[i][j] = exp(crf.edgePot[i][j]);

	UNPROTECT(crf.nEdges + 3);

	return(_crf);
}
