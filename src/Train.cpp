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

	double *nodeFea = NUMERIC_POINTER(AS_NUMERIC(_nodeFea));
	if (!ISNA(nodeFea[0]))
	{
		for (int i = 0; i < nNodes; i++)
		{
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
		}
	}

	double *edgeFea = NUMERIC_POINTER(AS_NUMERIC(_edgeFea));
	if (!ISNA(edgeFea[0]))
	{
		for (int i = 0; i < nEdges; i++)
		{
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
		}
	}

	if (Rf_isNewList(_nodeExt))
	{
		for (int i = 0; i < nPar; i++)
		{
			_temp = AS_NUMERIC(VECTOR_ELT(_nodeExt, i));
			double *nodeExt = NUMERIC_POINTER(_temp);
			if (!ISNA(nodeExt[0]))
			{
				for (int j = 0; j < nNodes; j++)
				{
					for (int k = 0; k < nStates[j]; k++)
					{
						nodePot[j + nNodes * k] += nodeExt[j + nNodes * k] * par[i];
					}
				}
			}
		}
	}

	if (Rf_isNewList(_edgeExt))
	{
		for (int i = 0; i < nPar; i++)
		{
			_temp = AS_LIST(VECTOR_ELT(_edgeExt, i));
			if (Rf_isNewList(_temp))
			{
				for (int j = 0; j < nEdges; j++)
				{
					double *edgeExt = NUMERIC_POINTER(AS_NUMERIC(VECTOR_ELT(_temp, j)));
					if (!ISNA(edgeExt[0]))
					{
						for (int k = 0; k < nEdgeStates[j]; k++)
						{
							edgePot[j][k] += edgeExt[k] * par[i];
						}
					}
				}
			}
		}
	}

	for (int i = 0; i < nNodes * maxState; i++)
		nodePot[i] = exp(nodePot[i]);
	for (int i = 0; i < nEdges; i++)
		for (int j = 0; j < nEdgeStates[i]; j++)
			edgePot[i][j] = exp(edgePot[i][j]);

	UNPROTECT(nEdges + 3);
}

void CRF::Update_Pot()
{
	int nPar = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.par")))[0];

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
	{
		for (int k = 0; k < nStates[i]; k++)
		{
			int p = nodePar[i + nNodes * k] - 1;
			if (p >= 0 && p < nPar)
				nodePot[i + nNodes * k] += par[p];
		}
	}

	for (int i = 0; i < nEdges; i++)
	{
		for (int k = 0; k < nEdgeStates[i]; k++)
		{
			int p = edgePar[i][k] - 1;
			if (p >= 0 && p < nPar)
				edgePot[i][k] += par[p];
		}
	}

	for (int i = 0; i < nNodes * maxState; i++)
		nodePot[i] = exp(nodePot[i]);
	for (int i = 0; i < nEdges; i++)
		for (int j = 0; j < nEdgeStates[i]; j++)
			edgePot[i][j] = exp(edgePot[i][j]);

	UNPROTECT(nEdges + 3);
}

SEXP MRF_Stat(SEXP _crf, SEXP _nInstances, SEXP _instances)
{
	CRF crf(_crf);

	int nInstances = INTEGER_POINTER(AS_INTEGER(_nInstances))[0];
	int nPar = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.par")))[0];

	double *instances = NUMERIC_POINTER(AS_NUMERIC(_instances));

	SEXP _stat;
	double *stat;
	PROTECT(_stat = NEW_NUMERIC(nPar));
	stat = NUMERIC_POINTER(_stat);
	SetValues(_stat, stat, 0.0);

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
	{
		for (int i = 0; i < crf.nNodes; i++)
		{
			int k = instances[n + nInstances * i] - 1;
			int p = nodePar[i + crf.nNodes * k] - 1;
			if (p >= 0 && p < nPar)
				stat[p]++;
		}

		for (int i = 0; i < crf.nEdges; i++)
		{
			int s1 = instances[n + nInstances * crf.EdgesBegin(i)] - 1;
			int s2 = instances[n + nInstances * crf.EdgesEnd(i)] - 1;
			int k = s1 + crf.nStates[crf.EdgesBegin(i)] * s2;
			int p = edgePar[i][k] - 1;
			if (p >= 0 && p < nPar)
				stat[p]++;
		}
	}

	UNPROTECT(crf.nEdges + 3);

	return(_stat);
}

SEXP MRF_NLL(SEXP _crf, SEXP _par, SEXP _nInstances, SEXP _infer, SEXP _env)
{
	CRF crf(_crf);

	int nInstances = INTEGER_POINTER(AS_INTEGER(_nInstances))[0];
	int nPar = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.par")))[0];

	double *par = NUMERIC_POINTER(AS_NUMERIC(_par));

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

	SEXP _crfPar, _parStat;
	double *crfPar, *parStat;
	PROTECT(_crfPar = AS_NUMERIC(GetListElement(_crf, "par")));
	PROTECT(_parStat = AS_NUMERIC(GetListElement(_crf, "par.stat")));
	crfPar = NUMERIC_POINTER(_crfPar);
	parStat = NUMERIC_POINTER(_parStat);
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

	crf.Update_Pot();

	SEXP _belief, _nodeBel, _edgeBel;
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

	*nll = NUMERIC_POINTER(AS_NUMERIC(GetListElement(_belief, "logZ")))[0] * nInstances;
	for (int i = 0; i < nPar; i++)
	{
		*nll -= par[i] * parStat[i];
		gradient[i] = -parStat[i];
	}

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

	for (int i = 0; i < crf.nEdges; i++)
	{
		for (int k = 0; k < crf.nEdgeStates[i]; k++)
		{
			int p = edgePar[i][k] - 1;
			if (p >= 0 && p < nPar)
			{
				gradient[p] += edgeBel[i][k] * nInstances;
			}
		}
	}

	UNPROTECT(crf.nEdges * 2 + 9);

	return(_nll);
}

SEXP CRF_NLL(SEXP _crf, SEXP _par, SEXP _nInstances, SEXP _instances, SEXP _nodeFea, SEXP _edgeFea, SEXP _nodeExt, SEXP _edgeExt, SEXP _infer, SEXP _env)
{
	CRF crf(_crf);

	int nInstances = INTEGER_POINTER(AS_INTEGER(_nInstances))[0];
	int nPar = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.par")))[0];
	int nNodeFea = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.nf")))[0];
	int nEdgeFea = INTEGER_POINTER(AS_INTEGER(GetListElement(_crf, "n.ef")))[0];

	double *par = NUMERIC_POINTER(AS_NUMERIC(_par));
	double *instances = NUMERIC_POINTER(AS_NUMERIC(_instances));

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

		SEXP _belief, _nodeBel, _edgeBel;
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

		double *nodeFea = NUMERIC_POINTER(AS_NUMERIC(_nodeFeaN));
		if (!ISNA(nodeFea[0]))
		{
			for (int i = 0; i < crf.nNodes; i++)
			{
				int s = y[i];
				for (int j = 0; j < nNodeFea; j++)
				{
					double f = nodeFea[j + nNodeFea * i];
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
		}

		double *edgeFea = NUMERIC_POINTER(AS_NUMERIC(_edgeFeaN));
		if (!ISNA(edgeFea[0]))
		{
			for (int i = 0; i < crf.nEdges; i++)
			{
				int s = y[crf.EdgesBegin(i)] + crf.nStates[crf.EdgesBegin(i)] * y[crf.EdgesEnd(i)];
				for (int j = 0; j < nEdgeFea; j++)
				{
					double f = edgeFea[j + nEdgeFea * i];
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
		}

		if (Rf_isNewList(_nodeExtN))
		{
			for (int i = 0; i < nPar; i++)
			{
				_temp = AS_NUMERIC(VECTOR_ELT(_nodeExtN, i));
				double *nodeExt = NUMERIC_POINTER(_temp);
				if (!ISNA(nodeExt[0]))
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
			}
		}

		if (Rf_isNewList(_edgeExt))
		{
			for (int i = 0; i < nPar; i++)
			{
				_temp = AS_LIST(VECTOR_ELT(_edgeExtN, i));
				if (Rf_isNewList(_temp))
				{
					for (int j = 0; j < crf.nEdges; j++)
					{
						double *edgeExt = NUMERIC_POINTER(AS_NUMERIC(VECTOR_ELT(_temp, j)));
						if (!ISNA(edgeExt[0]))
						{
							int s = y[crf.EdgesBegin(i)] + crf.nStates[crf.EdgesBegin(i)] * y[crf.EdgesEnd(i)];
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
					}
				}
			}
		}

		UNPROTECT(crf.nEdges + 3);
	}

	UNPROTECT(crf.nEdges + 5);

	return(_nll);
}
