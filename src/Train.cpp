#include "CRF.h"

SEXP Update_Pot(SEXP _crf, SEXP _nodeFea, SEXP _edgeFea)
{
	CRF crf(_crf);
	crf.Update_Pot(NUMERIC_POINTER(AS_NUMERIC(_nodeFea)), NUMERIC_POINTER(AS_NUMERIC(_edgeFea)));
	return (_crf);
}

void CRF::Update_Pot(double *nodeFea, double *edgeFea)
{
	int nPar = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.par")))[0];
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
	edgePar = (int **) R_alloc(nEdges, sizeof(int *));
	for (int i = 0; i < nEdges; i++)
	{
		PROTECT(_temp = AS_INTEGER(VECTOR_ELT(_edgePar, i)));
		edgePar[i] = INTEGER_POINTER(_temp);
	}

	for (int i = 0; i < nNodes * maxState; i++)
		nodePot[i] = 0;
	for (int i = 0; i < nEdges; i++)
		for (int j = 0; j < nEdgeStates[i]; j++)
			edgePot[i][j] = 0;

	for (int i = 0; i < nNodes; i++)
		for (int j = 0; j < nNodeFea; j++)
		{
			double f = nodeFea[j + nNodeFea * i];
			if (f != 0)
				for (int k = 0; k < nStates[i]; k++)
				{
					int p = nodePar[i + nNodes * (k + maxState * j)] - 1;
					if (p >= 0 && p < nPar)
						nodePot[i + nNodes * k] += f * par[p];
				}
		}

	for (int i = 0; i < nEdges; i++)
		for (int j = 0; j < nEdgeFea; j++)
		{
			double f = edgeFea[j + nEdgeFea * i];
			if (f != 0)
				for (int k = 0; k < nEdgeStates[i]; k++)
				{
					int p = edgePar[i][k + nEdgeStates[i] * j] - 1;
					if (p >= 0 && p < nPar)
						edgePot[i][k] += f * par[p];
				}
		}

	for (int i = 0; i < nNodes * maxState; i++)
		nodePot[i] = exp(nodePot[i]);
	for (int i = 0; i < nEdges; i++)
		for (int j = 0; j < nEdgeStates[i]; j++)
			edgePot[i][j] = exp(edgePot[i][j]);

	UNPROTECT(nEdges + 3);
}

SEXP Update_ParStat(SEXP _crf, SEXP _nInstances, SEXP _instances, SEXP _nodeFea, SEXP _edgeFea)
{
	CRF crf(_crf);

	int nInstances = INTEGER_POINTER(AS_INTEGER(_nInstances))[0];
	int nPar = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.par")))[0];
	int nNodeFea = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.nf")))[0];
	int nEdgeFea = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.ef")))[0];

	double *instances = NUMERIC_POINTER(AS_NUMERIC(_instances));
	double *nodeFea = NUMERIC_POINTER(AS_NUMERIC(_nodeFea));
	double *edgeFea = NUMERIC_POINTER(AS_NUMERIC(_edgeFea));

	SEXP _wPar;
	double *wPar;
	PROTECT(_wPar = AS_NUMERIC(GetListElement(_crf, "w.par")));
	wPar = NUMERIC_POINTER(_wPar);
	SetValues(_wPar, wPar, 0.0);

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
						wPar[p] += f;
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
						wPar[p] += f;
				}
			}
		}

	UNPROTECT(crf.nEdges + 3);

	return(_crf);
}

SEXP CRF_NLL(SEXP _crf, SEXP _par, SEXP _nInstances, SEXP _instances, SEXP _nodeFea, SEXP _edgeFea, SEXP _infer, SEXP _env)
{
	CRF crf(_crf);

	int nInstances = INTEGER_POINTER(AS_INTEGER(_nInstances))[0];
	int nPar = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.par")))[0];
	int nNodeFea = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.nf")))[0];
	int nEdgeFea = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.ef")))[0];

	double *par = NUMERIC_POINTER(AS_NUMERIC(_par));
	double *instances = NUMERIC_POINTER(AS_NUMERIC(_instances));
	double *nodeFea = NUMERIC_POINTER(AS_NUMERIC(_nodeFea));
	double *edgeFea = NUMERIC_POINTER(AS_NUMERIC(_edgeFea));

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

	SEXP _crfPar;
	double *crfPar;
	PROTECT(_crfPar = AS_NUMERIC(GetListElement(_crf, "par")));
	crfPar = NUMERIC_POINTER(_crfPar);
	for (int i = 0; i < nPar; i++)
		crfPar[i] = par[i];

	SEXP _nll, _gradient;
	double *nll, *gradient;
	PROTECT(_nll = AS_NUMERIC(GetListElement(_crf, "nll")));
	PROTECT(_gradient = AS_NUMERIC(GetListElement(_crf, "gradient")));
	nll = NUMERIC_POINTER(_nll);
	gradient = NUMERIC_POINTER(_gradient);
	*nll = 0.0;
	SetValues(_gradient, gradient, 0.0);

	int *y = (int *) R_allocVector<int>(crf.nNodes);

	for (int n = 0; n < nInstances; n++)
	{
		crf.Update_Pot(nodeFea + nNodeFea*crf.nNodes*n, edgeFea + nEdgeFea*crf.nEdges*n);

		for (int i = 0; i < crf.nNodes; i++)
			y[i] = instances[n + nInstances * i] - 1;

		SEXP _belief, _nodeBel, _edgeBel, _temp;
		double *nodeBel, **edgeBel;
		PROTECT(_belief = eval(_infer, _env));
		PROTECT(_nodeBel = AS_NUMERIC(GetListElement(_belief, "node.bel")));
		PROTECT(_edgeBel = AS_LIST(GetListElement(_belief, "edge.bel")));
		nodeBel = NUMERIC_POINTER(_nodeBel);
		edgeBel = (double **) R_alloc(crf.nEdges, sizeof(double *));
		for (int i = 0; i < crf.nEdges; i++)
		{
			PROTECT(_temp = AS_NUMERIC(VECTOR_ELT(_edgeBel, i)));
			edgeBel[i] = NUMERIC_POINTER(_temp);
		}

		*nll += NUMERIC_POINTER(AS_NUMERIC(GetListElement(_belief, "logZ")))[0] - crf.Get_LogPotential(y);

		for (int i = 0; i < crf.nNodes; i++)
		{
			int s = y[i];
			for (int j = 0; j < nNodeFea; j++)
			{
				double f = nodeFea[j + nNodeFea * (i + crf.nNodes * n)];
				if (f != 0)
				{
					for (int k = 0; k < crf.nStates[i]; k++)
					{
						int p = nodePar[i + crf.nNodes * (k + crf.maxState * j)] - 1;
						if (p >= 0 && p < nPar)
						{
							if (k == s)
							{
								gradient[p] -= f;
							}
							gradient[p] += f * nodeBel[i + crf.nNodes * k];
						}
					}
				}
			}
		}

		for (int i = 0; i < crf.nEdges; i++)
		{
			int s = y[crf.EdgesBegin(i)] + crf.nStates[crf.EdgesBegin(i)] * y[crf.EdgesEnd(i)];
			for (int j = 0; j < nEdgeFea; j++)
			{
				double f = edgeFea[j + nEdgeFea * (i + crf.nEdges * n)];
				if (f != 0)
				{
					for (int k = 0; k < crf.nEdgeStates[i]; k++)
					{
						int p = edgePar[i][k + crf.nEdgeStates[i] * j] - 1;
						if (p >= 0 && p < nPar)
						{
							if (k == s)
							{
								gradient[p] -= f;
							}
							gradient[p] += f * edgeBel[i][k];
						}
					}
				}
			}
		}

		UNPROTECT(crf.nEdges + 3);
	}

	UNPROTECT(crf.nEdges + 5);

	return(_nll);
}
