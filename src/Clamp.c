#include "CRF.h"

void Clamp_Reset(double *nodePotNew, int *clamped, int *nodeMap, int nNodesClamped, int nNodes, int nEdges, int *edges, int *nStates, int maxState, int *nAdj, int **adjEdges, double *nodePot, double *edgePot)
{
	int e, n1, n2;
	double *p_nodePot, *p_edgePot, *p_nodePotNew;

	for (int i = 0; i < nNodes; i++)
	{
		if (nodeMap[i] > 0)
		{
			p_nodePot = nodePot + i;
			p_nodePotNew = nodePotNew + nodeMap[i] - 1;
			for (int j = 0; j < nStates[i]; j++)
			{
				p_nodePotNew[0] = p_nodePot[0];
				p_nodePot += nNodes;
				p_nodePotNew += nNodesClamped;
			}
		}
	}

	for (int i = 0; i < nNodes; i++)
	{
		if (clamped[i])
		{
			for (int j = 0; j < nAdj[i]; j++)
			{
				e = adjEdges[i][j] - 1;
				n1 = edges[e] - 1;
				n2 = edges[e + nEdges] - 1;
				if (n1 == i && clamped[n2] == 0)
				{
					p_edgePot = edgePot + clamped[i] - 1 + maxState * maxState * e;
					p_nodePotNew = nodePotNew + nodeMap[n2] - 1;
					for (int k = 0; k < nStates[n2]; k++)
					{
						p_nodePotNew[0] *= p_edgePot[0];
						p_edgePot += maxState;
						p_nodePotNew += nNodesClamped;
					}
				}
				else if (n2 == i && clamped[n1] == 0)
				{
					p_edgePot = edgePot + maxState * (clamped[i] - 1 + maxState * e);
					p_nodePotNew = nodePotNew + nodeMap[n1] - 1;
					for (int k = 0; k < nStates[n1]; k++)
					{
						p_nodePotNew[0] *= p_edgePot[k];
						p_nodePotNew += nNodesClamped;
					}
				}
			}
		}
	}
}

SEXP Clamp_NodePot(SEXP _crfClamped)
{

	SEXP _clamped, _nodeMap, _nNodesClamped, _maxStateClamped;
	PROTECT(_clamped = AS_INTEGER(getListElement(_crfClamped, "clamped")));
	PROTECT(_nodeMap = AS_INTEGER(getListElement(_crfClamped, "node.map")));
	PROTECT(_nNodesClamped = AS_INTEGER(getListElement(_crfClamped, "n.nodes")));
	PROTECT(_maxStateClamped = AS_INTEGER(getListElement(_crfClamped, "max.state")));
	int *clamped = INTEGER_POINTER(_clamped);
	int *nodeMap = INTEGER_POINTER(_nodeMap);
	int nNodesClamped = INTEGER_POINTER(_nNodesClamped)[0];
	int maxStateClamped = INTEGER_POINTER(_maxStateClamped)[0];

	SEXP _crf;
	PROTECT(_crf = AS_LIST(getListElement(_crfClamped, "original")));

	SEXP _nNodes, _nEdges, _edges, _nStates, _maxState;
	PROTECT(_nNodes = AS_INTEGER(getListElement(_crf, "n.nodes")));
	PROTECT(_nEdges = AS_INTEGER(getListElement(_crf, "n.edges")));
	PROTECT(_edges = AS_INTEGER(getListElement(_crf, "edges")));
	PROTECT(_nStates = AS_INTEGER(getListElement(_crf, "n.states")));
	PROTECT(_maxState = AS_INTEGER(getListElement(_crf, "max.state")));
	int nNodes = INTEGER_POINTER(_nNodes)[0];
	int nEdges = INTEGER_POINTER(_nEdges)[0];
	int *edges = INTEGER_POINTER(_edges);
	int *nStates = INTEGER_POINTER(_nStates);
	int maxState = INTEGER_POINTER(_maxState)[0];

	SEXP _nAdj, _adjEdges, _temp;
	PROTECT(_nAdj = AS_INTEGER(getListElement(_crf, "n.adj")));
	PROTECT(_adjEdges = AS_LIST(getListElement(_crf, "adj.edges")));
	int *nAdj = INTEGER_POINTER(_nAdj);
	int **adjEdges = (int **) R_alloc(nNodes, sizeof(int *));
	for (int i = 0; i < nNodes; i++)
	{
		PROTECT(_temp = AS_INTEGER(VECTOR_ELT(_adjEdges, i)));
		adjEdges[i] = INTEGER_POINTER(_temp);
	}

	SEXP _nodePot, _edgePot;
	PROTECT(_nodePot = AS_NUMERIC(getListElement(_crf, "node.pot")));
	PROTECT(_edgePot = AS_NUMERIC(getListElement(_crf, "edge.pot")));
	double *nodePot = NUMERIC_POINTER(_nodePot);
	double *edgePot = NUMERIC_POINTER(_edgePot);

	SEXP _nodePotNew;
	PROTECT(_nodePotNew = NEW_NUMERIC(nNodesClamped * maxStateClamped));
	setDim2(_nodePotNew, nNodesClamped, maxStateClamped);
	double *nodePotNew = NUMERIC_POINTER(_nodePotNew);
	setValues(_nodePotNew, nodePotNew, 0);

	Clamp_Reset(nodePotNew, clamped, nodeMap, nNodesClamped, nNodes, nEdges, edges, nStates, maxState, nAdj, adjEdges, nodePot, edgePot);

	UNPROTECT(15 + nNodes);

	return(_nodePotNew);
}
