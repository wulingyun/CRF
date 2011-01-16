#include "CRF.h"

void CRF::Clamp_Reset(int *clamped, int *nodeMap, int nNodesNew, double *nodePotNew)
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
				p_nodePotNew += nNodesNew;
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
						p_nodePotNew += nNodesNew;
					}
				}
				else if (n2 == i && clamped[n1] == 0)
				{
					p_edgePot = edgePot + maxState * (clamped[i] - 1 + maxState * e);
					p_nodePotNew = nodePotNew + nodeMap[n1] - 1;
					for (int k = 0; k < nStates[n1]; k++)
					{
						p_nodePotNew[0] *= p_edgePot[k];
						p_nodePotNew += nNodesNew;
					}
				}
			}
		}
	}
}

SEXP Clamp_NodePot(SEXP _crfClamped)
{
	SEXP _crf;
	PROTECT(_crf = AS_LIST(getListElement(_crfClamped, "original")));
	CRF crf(_crf);

	SEXP _clamped, _nodeMap, _nNodes, _maxState;
	PROTECT(_clamped = AS_INTEGER(getListElement(_crfClamped, "clamped")));
	PROTECT(_nodeMap = AS_INTEGER(getListElement(_crfClamped, "node.map")));
	PROTECT(_nNodes = AS_INTEGER(getListElement(_crfClamped, "n.nodes")));
	PROTECT(_maxState = AS_INTEGER(getListElement(_crfClamped, "max.state")));
	int *clamped = INTEGER_POINTER(_clamped);
	int *nodeMap = INTEGER_POINTER(_nodeMap);
	int nNodes = INTEGER_POINTER(_nNodes)[0];
	int maxState = INTEGER_POINTER(_maxState)[0];

	SEXP _nodePot;
	PROTECT(_nodePot = NEW_NUMERIC(nNodes * maxState));
	setDim2(_nodePot, nNodes, maxState);
	double *nodePot = NUMERIC_POINTER(_nodePot);
	setValues(_nodePot, nodePot, 0);

	crf.Clamp_Reset(clamped, nodeMap, nNodes, nodePot);

	UNPROTECT(6);

	return(_nodePot);
}
