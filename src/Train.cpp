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
	int nPar = INTEGER_POINTER(GetVarAsInteger(_crf, "n.par"))[0];

	SEXP _par;
	PROTECT(_par = GetVarAsNumeric(_crf, "par"));
	double *par = NUMERIC_POINTER(_par);

	for (int i = 0; i < nNodes * maxState; i++)
		nodePot[i] = 0;
	for (int i = 0; i < nEdges; i++)
		for (int j = 0; j < nEdgeStates[i]; j++)
			edgePot[i][j] = 0;

  if (!isNull(_nodeFea))
  {
    PROTECT(_nodeFea = AS_NUMERIC(_nodeFea));
  	double *nodeFea = NUMERIC_POINTER(_nodeFea);
  	if (!ISNAN(nodeFea[0]))
  	{
  		int nNodeFea = INTEGER_POINTER(GetVarAsInteger(_crf, "n.nf"))[0];
  		SEXP _nodePar;
  		PROTECT(_nodePar = GetVarAsInteger(_crf, "node.par"));
  		int *nodePar = INTEGER_POINTER(_nodePar);
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
  		UNPROTECT(1);
  	}
    UNPROTECT(1);
  }

  if (!isNull(_edgeFea))
  {
  	PROTECT(_edgeFea = AS_NUMERIC(_edgeFea));
  	double *edgeFea = NUMERIC_POINTER(_edgeFea);
  	if (!ISNAN(edgeFea[0]))
  	{
  		int nEdgeFea = INTEGER_POINTER(GetVarAsInteger(_crf, "n.ef"))[0];
  		SEXP _edgePar;
  		PROTECT(_edgePar = GetVar(_crf, "edge.par"));
  		for (int i = 0; i < nEdges; i++)
  		{
  			SEXP _edgeParI;
  			PROTECT(_edgeParI = AS_INTEGER(GetListElement(_edgePar, i)));
  			int *edgePar = INTEGER_POINTER(_edgeParI);
  			for (int j = 0; j < nEdgeFea; j++)
  			{
  				double f = edgeFea[j + nEdgeFea * i];
  				if (f != 0)
  					for (int k = 0; k < nEdgeStates[i]; k++)
  					{
  						int p = edgePar[k + nEdgeStates[i] * j] - 1;
  						if (p >= 0 && p < nPar)
  							edgePot[i][k] += f * par[p];
  					}
  			}
  			UNPROTECT(1);
  		}
  		UNPROTECT(1);
  	}
    UNPROTECT(1);
  }

	if (!isNull(_nodeExt) && isNewList(_nodeExt))
	{
		for (int i = 0; i < nPar; i++)
		{
			SEXP _nodeExtI = GetListElement(_nodeExt, i);
      if (!isNull(_nodeExtI))
      {
  			PROTECT(_nodeExtI = AS_NUMERIC(_nodeExtI));
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
	}

	if (!isNull(_edgeExt) && isNewList(_edgeExt))
	{
		for (int i = 0; i < nPar; i++)
		{
			SEXP _edgeExtI = GetListElement(_edgeExt, i);
			if (!isNull(_edgeExtI) && isNewList(_edgeExtI))
			{
				for (int j = 0; j < nEdges; j++)
				{
					SEXP _edgeExtII = GetListElement(_edgeExtI, j);
          if (!isNull(_edgeExtII))
          {
  					PROTECT(_edgeExtII = AS_NUMERIC(_edgeExtII));
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
	}

	Update_Pot_Finalize();

	UNPROTECT(1);
}

void CRF::Update_Pot()
{
	int nPar = INTEGER_POINTER(GetVarAsInteger(_crf, "n.par"))[0];

	SEXP _par;
	PROTECT(_par = GetVarAsNumeric(_crf, "par"));
	double *par = NUMERIC_POINTER(_par);

	for (int i = 0; i < nNodes * maxState; i++)
		nodePot[i] = 0;
	for (int i = 0; i < nEdges; i++)
		for (int j = 0; j < nEdgeStates[i]; j++)
			edgePot[i][j] = 0;

	SEXP _nodePar;
	PROTECT(_nodePar = GetVarAsInteger(_crf, "node.par"));
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

	SEXP _edgePar;
	PROTECT(_edgePar = GetVar(_crf, "edge.par"));
	for (int i = 0; i < nEdges; i++)
	{
		SEXP _edgeParI;
		PROTECT(_edgeParI = AS_INTEGER(GetListElement(_edgePar, i)));
		int *edgePar = INTEGER_POINTER(_edgeParI);
		for (int k = 0; k < nEdgeStates[i]; k++)
		{
			int p = edgePar[k] - 1;
			if (p >= 0 && p < nPar)
				edgePot[i][k] += par[p];
		}
		UNPROTECT(1);
	}

	Update_Pot_Finalize();

	UNPROTECT(3);
}

void CRF::Update_Pot_Finalize()
{
  double maxPot;
  for (int i = 0; i < nNodes; i++)
  {
    maxPot = 0;
    for (int j = 0; j < nStates[i]; j++)
      maxPot = max(maxPot, NodePot(i, j));
    for (int j = 0; j < nStates[i]; j++)
      NodePot(i, j) -= maxPot;
  }
  
  int n1, n2;
  for (int i = 0; i < nEdges; i++)
  {
    n1 = EdgesBegin(i);
    n2 = EdgesEnd(i);
    maxPot = 0;
    for (int j = 0; j < nStates[n2]; j++)
    {
      for (int k = 0; k < nStates[n1]; k++)
        maxPot = max(maxPot, EdgePot(i, k, j));
    }
    for (int j = 0; j < nStates[n2]; j++)
    {
      for (int k = 0; k < nStates[n1]; k++)
        EdgePot(i, k, j) -= maxPot;
    }
  }
  
  for (int i = 0; i < nNodes * maxState; i++)
    nodePot[i] = max(exp(nodePot[i]), 1e-8);
  for (int i = 0; i < nEdges; i++)
    for (int j = 0; j < nEdgeStates[i]; j++)
      edgePot[i][j] = max(exp(edgePot[i][j]), 1e-8);
}

SEXP MRF_Stat(SEXP _crf, SEXP _instances)
{
	CRF crf(_crf);

	int nInstances = INTEGER_POINTER(GET_DIM(_instances))[0];
	int nPar = INTEGER_POINTER(GetVarAsInteger(_crf, "n.par"))[0];

	PROTECT(_instances = AS_NUMERIC(_instances));
	double *instances = NUMERIC_POINTER(_instances);

	SEXP _nodePar;
	PROTECT(_nodePar = GetVarAsInteger(_crf, "node.par"));
	int *nodePar = INTEGER_POINTER(_nodePar);

	SEXP _edgePar;
	PROTECT(_edgePar = GetVar(_crf, "edge.par"));
	int **edgePar = (int **) R_alloc(crf.nEdges, sizeof(int *));
	SEXP _edgeParI, _temp;
	PROTECT(_edgeParI = NEW_LIST(crf.nEdges));
	for (int i = 0; i < crf.nEdges; i++)
	{
		SET_VECTOR_ELT(_edgeParI, i, _temp = AS_INTEGER(GetListElement(_edgePar, i)));
		edgePar[i] = INTEGER_POINTER(_temp);
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

	UNPROTECT(5);

	return(_stat);
}

SEXP MRF_NLL(SEXP _crf, SEXP _par, SEXP _instances, SEXP _infer, SEXP _env)
{
	CRF crf(_crf);

	int nInstances = INTEGER_POINTER(GET_DIM(_instances))[0];
	int nPar = INTEGER_POINTER(GetVarAsInteger(_crf, "n.par"))[0];

	PROTECT(_par = AS_NUMERIC(_par));
	double *par = NUMERIC_POINTER(_par);
	SEXP _crfPar;
	PROTECT(_crfPar = GetVarAsNumeric(_crf, "par"));
	double *crfPar = NUMERIC_POINTER(_crfPar);
	for (int i = 0; i < nPar; i++)
		crfPar[i] = par[i];

	SEXP _parStat;
	PROTECT(_parStat = GetVarAsNumeric(_crf, "par.stat"));
	double *parStat = NUMERIC_POINTER(_parStat);

	SEXP _nll;
	PROTECT(_nll = GetVar(_crf, "nll"));
	double *nll = NUMERIC_POINTER(_nll);
	*nll = 0.0;

	SEXP _gradient;
	PROTECT(_gradient = GetVarAsNumeric(_crf, "gradient"));
	double *gradient = NUMERIC_POINTER(_gradient);
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
	PROTECT(_nodePar = GetVarAsInteger(_crf, "node.par"));
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

	SEXP _edgePar, _edgeBel;
	PROTECT(_edgePar = GetVar(_crf, "edge.par"));
	PROTECT(_edgeBel = GetListElement(_belief, "edge.bel"));
	SEXP _edgeParI, _edgeBelI, _temp;
	PROTECT(_edgeParI = NEW_LIST(crf.nEdges));
	PROTECT(_edgeBelI = NEW_LIST(crf.nEdges));
	for (int i = 0; i < crf.nEdges; i++)
	{
		SET_VECTOR_ELT(_edgeParI, i, _temp = AS_INTEGER(GetListElement(_edgePar, i)));
	  int *edgePar = INTEGER_POINTER(_temp);
	  SET_VECTOR_ELT(_edgeBelI, i, _temp = AS_NUMERIC(GetListElement(_edgeBel, i)));
		double *edgeBel = NUMERIC_POINTER(_temp);
		for (int k = 0; k < crf.nEdgeStates[i]; k++)
		{
			int p = edgePar[k] - 1;
			if (p >= 0 && p < nPar)
			{
				gradient[p] += edgeBel[k] * nInstances;
			}
		}
	}

	UNPROTECT(12);

	return(_nll);
}

SEXP CRF_NLL(SEXP _crf, SEXP _par, SEXP _instances, SEXP _nodeFea, SEXP _edgeFea, SEXP _nodeExt, SEXP _edgeExt, SEXP _infer, SEXP _env)
{
	CRF crf(_crf);

	int nInstances = INTEGER_POINTER(GET_DIM(_instances))[0];
	int nPar = INTEGER_POINTER(GetVarAsInteger(_crf, "n.par"))[0];
	int nNodeFea = INTEGER_POINTER(GetVarAsInteger(_crf, "n.nf"))[0];
	int nEdgeFea = INTEGER_POINTER(GetVarAsInteger(_crf, "n.ef"))[0];

	PROTECT(_par = AS_NUMERIC(_par));
	double *par = NUMERIC_POINTER(_par);
	SEXP _crfPar;
	PROTECT(_crfPar = GetVarAsNumeric(_crf, "par"));
	double *crfPar = NUMERIC_POINTER(_crfPar);
	for (int i = 0; i < nPar; i++)
		crfPar[i] = par[i];

	PROTECT(_instances = AS_NUMERIC(_instances));
	double *instances = NUMERIC_POINTER(_instances);

	SEXP _nodePar;
	PROTECT(_nodePar = GetVarAsInteger(_crf, "node.par"));
	int *nodePar = INTEGER_POINTER(_nodePar);

	SEXP _edgePar;
	PROTECT(_edgePar = GetVar(_crf, "edge.par"));
	int **edgePar = (int **) R_alloc(crf.nEdges, sizeof(int *));
	SEXP _edgeParI, _temp;
	PROTECT(_edgeParI = NEW_LIST(crf.nEdges));
	for (int i = 0; i < crf.nEdges; i++)
	{
		SET_VECTOR_ELT(_edgeParI, i, _temp = AS_INTEGER(GetListElement(_edgePar, i)));
		edgePar[i] = INTEGER_POINTER(_temp);
	}

	SEXP _nll;
	PROTECT(_nll = GetVar(_crf, "nll"));
	double *nll = NUMERIC_POINTER(_nll);
	*nll = 0.0;

	SEXP _gradient;
	PROTECT(_gradient = GetVarAsNumeric(_crf, "gradient"));
	double *gradient = NUMERIC_POINTER(_gradient);
	for (int i = 0; i < nPar; i++)
		gradient[i] = 0.0;

	int *y = (int *) R_allocVector<int>(crf.nNodes);

	SEXP _nodeFeaN = _nodeFea;
	SEXP _edgeFeaN = _edgeFea;
	SEXP _nodeExtN = _nodeExt;
	SEXP _edgeExtN = _edgeExt;
	for (int n = 0; n < nInstances; n++)
	{
		if (!isNull(_nodeFea) && isNewList(_nodeFea)) _nodeFeaN = GetListElement(_nodeFea, n);
		if (!isNull(_edgeFea) && isNewList(_edgeFea)) _edgeFeaN = GetListElement(_edgeFea, n);
		if (!isNull(_nodeExt) && isNewList(_nodeExt)) _nodeExtN = GetListElement(_nodeExt, n);
		if (!isNull(_edgeExt) && isNewList(_edgeExt)) _edgeExtN = GetListElement(_edgeExt, n);

		crf.Update_Pot(_nodeFeaN, _edgeFeaN, _nodeExtN, _edgeExtN);

		for (int i = 0; i < crf.nNodes; i++)
			y[i] = instances[n + nInstances * i] - 1;

		SEXP _belief;
		PROTECT(_belief = eval(_infer, _env));

		SEXP _nodeBel;
		PROTECT(_nodeBel = AS_NUMERIC(GetListElement(_belief, "node.bel")));
		double *nodeBel = NUMERIC_POINTER(_nodeBel);

		SEXP _edgeBel;
		PROTECT(_edgeBel = GetListElement(_belief, "edge.bel"));
		double **edgeBel = (double **) R_alloc(crf.nEdges, sizeof(double *));
		SEXP _edgeBelI, _temp;
		PROTECT(_edgeBelI = NEW_LIST(crf.nEdges));
		for (int i = 0; i < crf.nEdges; i++)
		{
			SET_VECTOR_ELT(_edgeBelI, i, _temp = AS_NUMERIC(GetListElement(_edgeBel, i)));
			edgeBel[i] = NUMERIC_POINTER(_temp);
		}

		*nll += NUMERIC_POINTER(AS_NUMERIC(GetListElement(_belief, "logZ")))[0] - crf.Get_LogPotential(y);

    if (!isNull(_nodeFeaN))
    {
  		PROTECT(_nodeFeaN = AS_NUMERIC(_nodeFeaN));
  		double *nodeFea = NUMERIC_POINTER(_nodeFeaN);
  		if (!ISNAN(nodeFea[0]))
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
      UNPROTECT(1);
    }

    if (!isNull(_edgeFeaN))
    {
  		PROTECT(_edgeFeaN = AS_NUMERIC(_edgeFeaN));
  		double *edgeFea = NUMERIC_POINTER(_edgeFeaN);
  		if (!ISNAN(edgeFea[0]))
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
      UNPROTECT(1);
    }

		if (!isNull(_nodeExtN) && isNewList(_nodeExtN))
		{
			for (int i = 0; i < nPar; i++)
			{
				SEXP _nodeExtI = GetListElement(_nodeExtN, i);
        if (!isNull(_nodeExtI))
        {
  				PROTECT(_nodeExtI = AS_NUMERIC(_nodeExtI));
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
		}

		if (!isNull(_edgeExtN) && isNewList(_edgeExtN))
		{
			for (int i = 0; i < nPar; i++)
			{
				SEXP _edgeExtI = GetListElement(_edgeExtN, i);
				if (!isNull(_edgeExtI) && isNewList(_edgeExtI))
				{
					for (int j = 0; j < crf.nEdges; j++)
					{
						SEXP _edgeExtII = GetListElement(_edgeExtI, j);
            if (!isNull(_edgeExtII))
            {
  						PROTECT(_edgeExtII = AS_NUMERIC(_edgeExtII));
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
		}

		UNPROTECT(4);
	}

	UNPROTECT(7);

	return(_nll);
}
