#include "CRF.h"

SEXP MRF_Update(SEXP _crf)
{
	CRF crf(_crf);
	crf.Update_Pot();
	return (_crf);
}

SEXP CRF_Update(SEXP _crf, SEXP _nodeFea, SEXP _edgeFea, SEXP _nodeExt, SEXP _edgeExt)
{
	CRF crf(_crf);
	crf.Update_Pot(_nodeFea, _edgeFea, _nodeExt, _edgeExt);
	return (_crf);
}

void CRF::Update_Pot(SEXP _nodeFea, SEXP _edgeFea, SEXP _nodeExt, SEXP _edgeExt)
{
	int nPar = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.par")))[0];

	SEXP _par;
	PROTECT(_par = AS_NUMERIC(GetListElement(_crf, "par")));
	double *par = NUMERIC_POINTER(_par);

	for (int i = 0; i < nNodes * maxState; i++)
		nodePot[i] = 0;
	for (int i = 0; i < nEdges; i++)
		for (int j = 0; j < nEdgeStates[i]; j++)
			edgePot[i][j] = 0;

	PROTECT(_nodeFea = AS_NUMERIC(_nodeFea));
	double *nodeFea = NUMERIC_POINTER(_nodeFea);
	if (!ISNAN(nodeFea[0]))
	{
		int nNodeFea = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.nf")))[0];
		SEXP _nodePar = GetListElement(_crf, "node.par");
		for (int i = 0; i < nNodeFea; i++)
		{
			SEXP _nodeParI;
			PROTECT(_nodeParI = AS_INTEGER(VECTOR_ELT(_nodePar, i)));
			int *nodePar = INTEGER_POINTER(_nodeParI);
			for (int j = 0; j < nNodes; j++)
			{
				double f = nodeFea[i + nNodeFea * j];
				if (f != 0)
				{
					for (int k = 0; k < nStates[j]; k++)
					{
						int p = nodePar[j + nNodes * k] - 1;
						if (p >= 0 && p < nPar)
							nodePot[j + nNodes * k] += f * par[p];
					}
				}
			}
			UNPROTECT(1);
		}
	}

	PROTECT(_edgeFea = AS_NUMERIC(_edgeFea));
	double *edgeFea = NUMERIC_POINTER(_edgeFea);
	if (!ISNAN(edgeFea[0]))
	{
		int nEdgeFea = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.ef")))[0];
		SEXP _edgePar = GetListElement(_crf, "edge.par");
		for (int i = 0; i < nEdgeFea; i++)
		{
			SEXP _edgeParI = VECTOR_ELT(_edgePar, i);
			for (int j = 0; j < nEdges; j++)
			{
				double f = edgeFea[i + nEdgeFea * j];
				if (f != 0)
				{
					SEXP _edgeParII;
					PROTECT(_edgeParII = AS_INTEGER(VECTOR_ELT(_edgeParI, j)));
					int *edgePar = INTEGER_POINTER(_edgeParII);
					for (int k = 0; k < nEdgeStates[j]; k++)
					{
						int p = edgePar[k] - 1;
						if (p >= 0 && p < nPar)
							edgePot[j][k] += f * par[p];
					}
					UNPROTECT(1);
				}
			}
		}
	}

	if (Rf_isNewList(_nodeExt))
	{
		for (int i = 0; i < nPar; i++)
		{
			SEXP _nodeExtI;
			PROTECT(_nodeExtI = AS_NUMERIC(VECTOR_ELT(_nodeExt, i)));
			double *nodeExt = NUMERIC_POINTER(_nodeExtI);
			if (!ISNAN(nodeExt[0]))
			{
				for (int j = 0; j < nNodes; j++)
				{
					for (int k = 0; k < nStates[j]; k++)
					{
						nodePot[j + nNodes * k] += nodeExt[j + nNodes * k] * par[i];
					}
				}
			}
			UNPROTECT(1);
		}
	}

	if (Rf_isNewList(_edgeExt))
	{
		for (int i = 0; i < nPar; i++)
		{
			SEXP _edgeExtI = VECTOR_ELT(_edgeExt, i);
			if (Rf_isNewList(_edgeExtI))
			{
				for (int j = 0; j < nEdges; j++)
				{
					SEXP _edgeExtII;
					PROTECT(_edgeExtII = AS_NUMERIC(VECTOR_ELT(_edgeExtI, j)));
					double *edgeExt = NUMERIC_POINTER(_edgeExtII);
					if (!ISNAN(edgeExt[0]))
					{
						for (int k = 0; k < nEdgeStates[j]; k++)
						{
							edgePot[j][k] += edgeExt[k] * par[i];
						}
					}
					UNPROTECT(1);
				}
			}
		}
	}

	for (int i = 0; i < nNodes * maxState; i++)
		nodePot[i] = exp(nodePot[i]);
	for (int i = 0; i < nEdges; i++)
		for (int j = 0; j < nEdgeStates[i]; j++)
			edgePot[i][j] = exp(edgePot[i][j]);

	UNPROTECT(3);
}

void CRF::Update_Pot()
{
	int nPar = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.par")))[0];

	SEXP _par;
	PROTECT(_par = AS_NUMERIC(GetListElement(_crf, "par")));
	double *par = NUMERIC_POINTER(_par);

	for (int i = 0; i < nNodes * maxState; i++)
		nodePot[i] = 0;
	for (int i = 0; i < nEdges; i++)
		for (int j = 0; j < nEdgeStates[i]; j++)
			edgePot[i][j] = 0;

	SEXP _nodePar;
	PROTECT(_nodePar = AS_INTEGER(VECTOR_ELT(GetListElement(_crf, "node.par"), 0)));
	int *nodePar = INTEGER_POINTER(_nodePar);
	for (int i = 0; i < nNodes; i++)
	{
		for (int k = 0; k < nStates[i]; k++)
		{
			int p = nodePar[i + nNodes * k] - 1;
			if (p >= 0 && p < nPar)
				nodePot[i + nNodes * k] += par[p];
		}
	}

	SEXP _edgePar = VECTOR_ELT(GetListElement(_crf, "edge.par"), 0);
	for (int i = 0; i < nEdges; i++)
	{
		SEXP _edgeParI;
		PROTECT(_edgeParI = AS_INTEGER(VECTOR_ELT(_edgePar, i)));
		int *edgePar = INTEGER_POINTER(_edgeParI);
		for (int k = 0; k < nEdgeStates[i]; k++)
		{
			int p = edgePar[k] - 1;
			if (p >= 0 && p < nPar)
				edgePot[i][k] += par[p];
		}
		UNPROTECT(1);
	}

	for (int i = 0; i < nNodes * maxState; i++)
		nodePot[i] = exp(nodePot[i]);
	for (int i = 0; i < nEdges; i++)
		for (int j = 0; j < nEdgeStates[i]; j++)
			edgePot[i][j] = exp(edgePot[i][j]);

	UNPROTECT(2);
}

SEXP MRF_Stat(SEXP _crf, SEXP _instances)
{
	CRF crf(_crf);

	int nInstances = INTEGER_POINTER(GET_DIM(_instances))[0];
	int nPar = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.par")))[0];

	PROTECT(_instances = AS_NUMERIC(_instances));
	double *instances = NUMERIC_POINTER(_instances);

	SEXP _nodePar;
	PROTECT(_nodePar = AS_INTEGER(VECTOR_ELT(GetListElement(_crf, "node.par"), 0)));
	int *nodePar = INTEGER_POINTER(_nodePar);

	SEXP _edgePar = VECTOR_ELT(GetListElement(_crf, "edge.par"), 0);
	int **edgePar = (int **) R_alloc(crf.nEdges, sizeof(int *));
	for (int i = 0; i < crf.nEdges; i++)
	{
		SEXP _edgeParI;
		PROTECT(_edgeParI = AS_INTEGER(VECTOR_ELT(_edgePar, i)));
		edgePar[i] = INTEGER_POINTER(_edgeParI);
	}

	SEXP _stat;
	PROTECT(_stat = NEW_NUMERIC(nPar));
	double *stat = NUMERIC_POINTER(_stat);
	SetValues(_stat, stat, 0.0);

	int *y = (int *) R_allocVector<int>(crf.nNodes);

	for (int n = 0; n < nInstances; n++)
	{
		for (int i = 0; i < crf.nNodes; i++)
		{
			y[i] = instances[n + nInstances * i] - 1;
			int p = nodePar[i + crf.nNodes * y[i]] - 1;
			if (p >= 0 && p < nPar)
				stat[p]++;
		}

		for (int i = 0; i < crf.nEdges; i++)
		{
			int p = edgePar[i][y[crf.EdgesBegin(i)] + crf.nStates[crf.EdgesBegin(i)] * y[crf.EdgesEnd(i)]] - 1;
			if (p >= 0 && p < nPar)
				stat[p]++;
		}
	}

	UNPROTECT(crf.nEdges + 3);

	return(_stat);
}

SEXP MRF_NLL(SEXP _crf, SEXP _par, SEXP _instances, SEXP _infer, SEXP _env)
{
	CRF crf(_crf);

	int nInstances = INTEGER_POINTER(GET_DIM(_instances))[0];
	int nPar = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.par")))[0];

	PROTECT(_par = AS_NUMERIC(_par));
	double *par = NUMERIC_POINTER(_par);
	double *crfPar = NUMERIC_POINTER(GetListElement(_crf, "par"));
	for (int i = 0; i < nPar; i++)
		crfPar[i] = par[i];

	SEXP _parStat;
	PROTECT(_parStat = AS_NUMERIC(GetListElement(_crf, "par.stat")));
	double *parStat = NUMERIC_POINTER(_parStat);

	SEXP _nll = GetListElement(_crf, "nll");
	double *nll = NUMERIC_POINTER(_nll);
	*nll = 0.0;

	double *gradient = NUMERIC_POINTER(GetListElement(_crf, "gradient"));
	for (int i = 0; i < nPar; i++)
		gradient[i] = 0.0;

	crf.Update_Pot();

	SEXP _belief;
	PROTECT(_belief = eval(_infer, _env));

	*nll = NUMERIC_POINTER(AS_NUMERIC(GetListElement(_belief, "logZ")))[0] * nInstances;
	for (int i = 0; i < nPar; i++)
	{
		*nll -= par[i] * parStat[i];
		gradient[i] = -parStat[i];
	}

	SEXP _nodePar, _nodeBel;
	PROTECT(_nodePar = AS_INTEGER(VECTOR_ELT(GetListElement(_crf, "node.par"), 0)));
	PROTECT(_nodeBel = AS_NUMERIC(GetListElement(_belief, "node.bel")));
	int *nodePar = INTEGER_POINTER(_nodePar);
	double *nodeBel = NUMERIC_POINTER(_nodeBel);
	for (int i = 0; i < crf.nNodes; i++)
	{
		for (int k = 0; k < crf.nStates[i]; k++)
		{
			int p = nodePar[i + crf.nNodes * k] - 1;
			if (p >= 0 && p < nPar)
			{
				gradient[p] += nodeBel[i + crf.nNodes * k] * nInstances;
			}
		}
	}

	SEXP _edgePar = VECTOR_ELT(GetListElement(_crf, "edge.par"), 0);
	SEXP _edgeBel = GetListElement(_belief, "edge.bel");
	for (int i = 0; i < crf.nEdges; i++)
	{
		SEXP _edgeParN, _edgeBelN;
		PROTECT(_edgeParN = AS_INTEGER(VECTOR_ELT(_edgePar, i)));
		PROTECT(_edgeBelN = AS_NUMERIC(VECTOR_ELT(_edgeBel, i)));
		int *edgePar = INTEGER_POINTER(_edgeParN);
		double *edgeBel = NUMERIC_POINTER(_edgeBelN);
		for (int k = 0; k < crf.nEdgeStates[i]; k++)
		{
			int p = edgePar[k] - 1;
			if (p >= 0 && p < nPar)
			{
				gradient[p] += edgeBel[k] * nInstances;
			}
		}
		UNPROTECT(2);
	}

	UNPROTECT(5);

	return(_nll);
}

SEXP CRF_NLL(SEXP _crf, SEXP _par, SEXP _instances, SEXP _nodeFea, SEXP _edgeFea, SEXP _nodeExt, SEXP _edgeExt, SEXP _infer, SEXP _env)
{
	CRF crf(_crf);

	int nInstances = INTEGER_POINTER(GET_DIM(_instances))[0];
	int nPar = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.par")))[0];
	int nNodeFea = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.nf")))[0];
	int nEdgeFea = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.ef")))[0];

	PROTECT(_par = AS_NUMERIC(_par));
	double *par = NUMERIC_POINTER(_par);
	double *crfPar = NUMERIC_POINTER(GetListElement(_crf, "par"));
	for (int i = 0; i < nPar; i++)
		crfPar[i] = par[i];

	PROTECT(_instances = AS_NUMERIC(_instances));
	double *instances = NUMERIC_POINTER(_instances);

	SEXP _nodePar = GetListElement(_crf, "node.par");
	int **nodePar = (int **) R_allocVector<int *>(nNodeFea);
	for (int i = 0; i < nNodeFea; i++)
	{
		SEXP _nodeParI;
		PROTECT(_nodeParI = AS_INTEGER(VECTOR_ELT(_nodePar, i)));
		nodePar[i] = INTEGER_POINTER(_nodeParI);
	}

	SEXP _edgePar = GetListElement(_crf, "edge.par");
	int ***edgePar = (int ***) R_allocArray<int *>(nEdgeFea, crf.nEdges);
	for (int i = 0; i < nEdgeFea; i++)
	{
		SEXP _edgeParI = VECTOR_ELT(_edgePar, i);
		for (int j = 0; j < crf.nEdges; j++)
		{
			SEXP _edgeParII;
			PROTECT(_edgeParII = AS_INTEGER(VECTOR_ELT(_edgeParI, j)));
			edgePar[i][j] = INTEGER_POINTER(_edgeParII);
		}
	}

	SEXP _nll = GetListElement(_crf, "nll");
	double *nll = NUMERIC_POINTER(_nll);
	*nll = 0.0;

	double *gradient = NUMERIC_POINTER(GetListElement(_crf, "gradient"));
	for (int i = 0; i < nPar; i++)
		gradient[i] = 0.0;

	int *y = (int *) R_allocVector<int>(crf.nNodes);

	SEXP _nodeFeaN = _nodeFea;
	SEXP _edgeFeaN = _edgeFea;
	SEXP _nodeExtN = _nodeExt;
	SEXP _edgeExtN = _edgeExt;
	for (int n = 0; n < nInstances; n++)
	{
		if (Rf_isNewList(_nodeFea)) _nodeFeaN = VECTOR_ELT(_nodeFea, n);
		if (Rf_isNewList(_edgeFea)) _edgeFeaN = VECTOR_ELT(_edgeFea, n);
		if (Rf_isNewList(_nodeExt)) _nodeExtN = VECTOR_ELT(_nodeExt, n);
		if (Rf_isNewList(_edgeExt)) _edgeExtN = VECTOR_ELT(_edgeExt, n);

		crf.Update_Pot(_nodeFeaN, _edgeFeaN, _nodeExtN, _edgeExtN);

		for (int i = 0; i < crf.nNodes; i++)
			y[i] = instances[n + nInstances * i] - 1;

		SEXP _belief;
		PROTECT(_belief = eval(_infer, _env));

		SEXP _nodeBel;
		PROTECT(_nodeBel = AS_NUMERIC(GetListElement(_belief, "node.bel")));
		double *nodeBel = NUMERIC_POINTER(_nodeBel);

		SEXP _edgeBel = GetListElement(_belief, "edge.bel");
		double **edgeBel = (double **) R_alloc(crf.nEdges, sizeof(double *));
		for (int i = 0; i < crf.nEdges; i++)
		{
			SEXP _edgeBelI;
			PROTECT(_edgeBelI = AS_NUMERIC(VECTOR_ELT(_edgeBel, i)));
			edgeBel[i] = NUMERIC_POINTER(_edgeBelI);
		}

		*nll += NUMERIC_POINTER(AS_NUMERIC(GetListElement(_belief, "logZ")))[0] - crf.Get_LogPotential(y);

		PROTECT(_nodeFeaN = AS_NUMERIC(_nodeFeaN));
		double *nodeFea = NUMERIC_POINTER(_nodeFeaN);
		if (!ISNAN(nodeFea[0]))
		{
			for (int i = 0; i < nNodeFea; i++)
			{
				for (int j = 0; j < crf.nNodes; j++)
				{
					int s = y[j];
					double f = nodeFea[i + nNodeFea * j];
					if (f != 0)
					{
						for (int k = 0; k < crf.nStates[j]; k++)
						{
							int p = nodePar[i][j + crf.nNodes * k] - 1;
							if (p >= 0 && p < nPar)
							{
								if (k == s)
								{
									gradient[p] -= f;
								}
								gradient[p] += f * nodeBel[j + crf.nNodes * k];
							}
						}
					}
				}
			}
		}

		PROTECT(_edgeFeaN = AS_NUMERIC(_edgeFeaN));
		double *edgeFea = NUMERIC_POINTER(_edgeFeaN);
		if (!ISNAN(edgeFea[0]))
		{
			for (int i = 0; i < nEdgeFea; i++)
			{
				for (int j = 0; j < crf.nEdges; j++)
				{
					int s = y[crf.EdgesBegin(j)] + crf.nStates[crf.EdgesBegin(j)] * y[crf.EdgesEnd(j)];
					double f = edgeFea[i + nEdgeFea * j];
					if (f != 0)
					{
						for (int k = 0; k < crf.nEdgeStates[j]; k++)
						{
							int p = edgePar[i][j][k] - 1;
							if (p >= 0 && p < nPar)
							{
								if (k == s)
								{
									gradient[p] -= f;
								}
								gradient[p] += f * edgeBel[j][k];
							}
						}
					}
				}
			}
		}

		if (Rf_isNewList(_nodeExtN))
		{
			for (int i = 0; i < nPar; i++)
			{
				SEXP _nodeExtI;
				PROTECT(_nodeExtI = AS_NUMERIC(VECTOR_ELT(_nodeExtN, i)));
				double *nodeExt = NUMERIC_POINTER(_nodeExtI);
				if (!ISNAN(nodeExt[0]))
				{
					for (int j = 0; j < crf.nNodes; j++)
					{
						int s = y[j];
						for (int k = 0; k < crf.nStates[j]; k++)
						{
							double f = nodeExt[j + crf.nNodes * k];
							if (k == s)
							{
								gradient[i] -= f;
							}
							gradient[i] += f * nodeBel[j + crf.nNodes * k];
						}
					}
				}
				UNPROTECT(1);
			}
		}

		if (Rf_isNewList(_edgeExt))
		{
			for (int i = 0; i < nPar; i++)
			{
				SEXP _edgeExtI = VECTOR_ELT(_edgeExtN, i);
				if (Rf_isNewList(_edgeExtI))
				{
					for (int j = 0; j < crf.nEdges; j++)
					{
						SEXP _edgeExtII;
						PROTECT(_edgeExtII = AS_NUMERIC(VECTOR_ELT(_edgeExtI, j)));
						double *edgeExt = NUMERIC_POINTER(_edgeExtII);
						if (!ISNAN(edgeExt[0]))
						{
							int s = y[crf.EdgesBegin(j)] + crf.nStates[crf.EdgesBegin(j)] * y[crf.EdgesEnd(j)];
							for (int k = 0; k < crf.nEdgeStates[j]; k++)
							{
								double f = edgeExt[k];
								if (k == s)
								{
									gradient[i] -= f;
								}
								gradient[i] += f * edgeBel[j][k];
							}
						}
						UNPROTECT(1);
					}
				}
			}
		}

		UNPROTECT(crf.nEdges + 4);
	}

	UNPROTECT(crf.nEdges * nEdgeFea + nNodeFea + 2);

	return(_nll);
}
