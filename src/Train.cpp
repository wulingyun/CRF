#include "CRF.h"

SEXP Update_Pot(SEXP _crf, SEXP _nodeFea, SEXP _edgeFea)
{
	CRF crf(_crf);

	int nPar = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.par")))[0];
	int nNodeFea = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.nf")))[0];
	int nEdgeFea = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.ef")))[0];

	double *nodeFea = NUMERIC_POINTER(AS_NUMERIC(_nodeFea));
	double *edgeFea = NUMERIC_POINTER(AS_NUMERIC(_edgeFea));

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
			double f = nodeFea[j + nNodeFea * i];
			if (f != 0)
				for (int k = 0; k < crf.nStates[i]; k++)
				{
					int p = nodePar[i + crf.nNodes * (k + crf.maxState * j)] - 1;
					if (p >= 0 && p < nPar)
						crf.nodePot[i + crf.nNodes * k] += f * par[p];
				}
		}

	for (int i = 0; i < crf.nEdges; i++)
		for (int j = 0; j < nEdgeFea; j++)
		{
			double f = edgeFea[j + nEdgeFea * i];
			if (f != 0)
				for (int k = 0; k < crf.nEdgeStates[i]; k++)
				{
					int p = edgePar[i][k + crf.nEdgeStates[i] * j] - 1;
					if (p >= 0 && p < nPar)
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

SEXP Get_SuffStat(SEXP _crf, SEXP _nInstances, SEXP _instances, SEXP _nodeFea, SEXP _edgeFea)
{
	CRF crf(_crf);

	int nInstances = INTEGER_POINTER(AS_INTEGER(_nInstances))[0];
	int nPar = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.par")))[0];
	int nNodeFea = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.nf")))[0];
	int nEdgeFea = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.ef")))[0];

	double *instances = NUMERIC_POINTER(AS_NUMERIC(_instances));
	double *nodeFea = NUMERIC_POINTER(AS_NUMERIC(_nodeFea));
	double *edgeFea = NUMERIC_POINTER(AS_NUMERIC(_edgeFea));

	SEXP _suffstat;
	double *suffstat;
	PROTECT(_suffstat = NEW_NUMERIC(nPar));
	suffstat = NUMERIC_POINTER(_suffstat);
	SetValues(_suffstat, suffstat, 0.0);

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

	for (int n = 0; n < nInstances; n++)
		for (int i = 0; i < crf.nNodes; i++)
		{
			int k = instances[n + nInstances * i] - 1;
			for (int j = 0; j < nNodeFea; j++)
			{
				double f = nodeFea[j + nNodeFea * (i + crf.nNodes * n)];
				if (f != 0)
				{
					int p = nodePar[i + crf.nNodes * (k + crf.maxState * j)] - 1;
					if (p >= 0 && p < nPar)
						suffstat[p] += f;
				}
			}
		}

	for (int n = 0; n < nInstances; n++)
		for (int i = 0; i < crf.nEdges; i++)
		{
			int s1 = instances[n + nInstances * crf.EdgesBegin(i)] - 1;
			int s2 = instances[n + nInstances * crf.EdgesEnd(i)] - 1;
			for (int j = 0; j < nEdgeFea; j++)
			{
				double f = edgeFea[j + nEdgeFea * (i + crf.nEdges * n)];
				if (f != 0)
				{
					int k = s1 + crf.nStates[crf.EdgesBegin(i)] * s2;
					int p = edgePar[i][k + crf.nEdgeStates[i] * j] - 1;
					if (p >= 0 && p < nPar)
						suffstat[p] += f;
				}
			}
		}

	UNPROTECT(crf.nEdges + 3);

	return(_suffstat);
}
