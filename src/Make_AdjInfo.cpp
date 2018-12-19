#include "CRF.h"

SEXP Make_AdjInfo(SEXP _crf)
{
	SEXP _nNodes, _nEdges, _edges;
	int nNodes, nEdges, *edges;

	PROTECT(_nNodes = GetVarAsInteger(_crf, "n.nodes"));
	PROTECT(_nEdges = GetVarAsInteger(_crf, "n.edges"));
	PROTECT(_edges = GetVarAsInteger(_crf, "edges"));
	nNodes = INTEGER_POINTER(_nNodes)[0];
	nEdges = INTEGER_POINTER(_nEdges)[0];
	edges = INTEGER_POINTER(_edges);

	int *nAdjTemp, **adjNodesTemp, **adjEdgesTemp;
	nAdjTemp = (int *) R_allocVector<int>(nNodes);
	adjNodesTemp = (int **) R_allocArray<int>(nNodes, nNodes);
	adjEdgesTemp = (int **) R_allocArray<int>(nNodes, nNodes);

	int n1, n2;
	for (int i = 0; i < nNodes; i++)
		nAdjTemp[i] = 0;
	for (int i = 0; i < nEdges; i++)
	{
		n1 = edges[i] - 1;
		n2 = edges[i + nEdges] - 1;
		adjNodesTemp[n1][nAdjTemp[n1]] = n2;
		adjNodesTemp[n2][nAdjTemp[n2]] = n1;
		adjEdgesTemp[n1][nAdjTemp[n1]] = i;
		adjEdgesTemp[n2][nAdjTemp[n2]] = i;
		nAdjTemp[n1]++;
		nAdjTemp[n2]++;
	}
	for (int i = 0; i < nNodes; i++)
	{
		R_isort(adjNodesTemp[i], nAdjTemp[i]);
		R_isort(adjEdgesTemp[i], nAdjTemp[i]);
	}


	SEXP _nAdj, _adjNodes, _adjEdges, _temp;
	int *nAdj, *temp1, *temp2;
	PROTECT(_nAdj = NEW_INTEGER(nNodes));
	PROTECT(_adjNodes = NEW_LIST(nNodes));
	PROTECT(_adjEdges = NEW_LIST(nNodes));
	nAdj = INTEGER_POINTER(_nAdj);
	for (int i = 0; i < nNodes; i++)
	{
		nAdj[i] = nAdjTemp[i];
	  SET_VECTOR_ELT(_adjNodes, i, _temp = NEW_INTEGER(nAdj[i]));
	  temp1 = INTEGER_POINTER(_temp);
	  SET_VECTOR_ELT(_adjEdges, i, _temp = NEW_INTEGER(nAdj[i]));
		temp2 = INTEGER_POINTER(_temp);
		for (int j = 0; j < nAdj[i]; j++)
		{
			temp1[j] = adjNodesTemp[i][j] + 1;
			temp2[j] = adjEdgesTemp[i][j] + 1;
		}
	}

	SetVar(_crf, "n.adj", _nAdj);
	SetVar(_crf, "adj.nodes", _adjNodes);
	SetVar(_crf, "adj.edges", _adjEdges);

	UNPROTECT(6);

	return(_crf);
}
