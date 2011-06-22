#include "CRF.h"

SEXP Clamp_NodePot(SEXP _crf)
{
	CRFclamped crf(_crf);
	return(crf._nodePot);
}

CRFclamped::CRFclamped(SEXP _crf)
: CRF(_crf)
{
	PROTECT(_original = AS_LIST(getListElement(_crf, "original")));
	original.Set_Data(_original);

	PROTECT(_nodeId = AS_INTEGER(getListElement(_crf, "node.id")));
	PROTECT(_nodeMap = AS_INTEGER(getListElement(_crf, "node.map")));
	PROTECT(_edgeId = AS_INTEGER(getListElement(_crf, "edge.id")));
	PROTECT(_edgeMap = AS_INTEGER(getListElement(_crf, "edge.map")));
	nodeId = INTEGER_POINTER(_nodeId);
	nodeMap = INTEGER_POINTER(_nodeMap);
	edgeId = INTEGER_POINTER(_edgeId);
	edgeMap = INTEGER_POINTER(_edgeMap);

	SEXP _clamped0;
	int *clamped0;
	PROTECT(_clamped0 = AS_INTEGER(getListElement(_crf, "clamped")));
	clamped0 = INTEGER_POINTER(_clamped0);

	PROTECT(_clamped = NEW_INTEGER(original.nNodes));
	clamped = INTEGER_POINTER(_clamped);
	for (int i = 0; i < original.nNodes; i++)
		clamped[i] = clamped0[i];

	PROTECT(_nodePot = NEW_NUMERIC(nNodes * maxState));
	setDim2(_nodePot, nNodes, maxState);
	nodePot = NUMERIC_POINTER(_nodePot);
	setValues(_nodePot, nodePot, 0);

	Reset_NodePot();

	numProtect += 8;
}

void CRFclamped::Reset_NodePot()
{
	int e, n1, n2;
	double *p_nodePot, *p_edgePot, *p_nodePotNew;

	for (int i = 0; i < original.nNodes; i++)
	{
		if (nodeMap[i] > 0)
		{
			p_nodePot = original.nodePot + i;
			p_nodePotNew = nodePot + nodeMap[i] - 1;
			for (int j = 0; j < original.nStates[i]; j++)
			{
				p_nodePotNew[0] = p_nodePot[0];
				p_nodePot += original.nNodes;
				p_nodePotNew += nNodes;
			}
		}
	}

	for (int i = 0; i < original.nNodes; i++)
	{
		if (clamped[i])
		{
			for (int j = 0; j < original.nAdj[i]; j++)
			{
				e = original.adjEdges[i][j] - 1;
				n1 = original.edges[e] - 1;
				n2 = original.edges[e + original.nEdges] - 1;
				if (n1 == i && clamped[n2] == 0)
				{
					p_edgePot = original.edgePot + clamped[i] - 1 + original.maxState * original.maxState * e;
					p_nodePotNew = nodePot + nodeMap[n2] - 1;
					for (int k = 0; k < original.nStates[n2]; k++)
					{
						p_nodePotNew[0] *= p_edgePot[0];
						p_edgePot += original.maxState;
						p_nodePotNew += nNodes;
					}
				}
				else if (n2 == i && clamped[n1] == 0)
				{
					p_edgePot = original.edgePot + original.maxState * (clamped[i] - 1 + original.maxState * e);
					p_nodePotNew = nodePot + nodeMap[n1] - 1;
					for (int k = 0; k < original.nStates[n1]; k++)
					{
						p_nodePotNew[0] *= p_edgePot[k];
						p_nodePotNew += nNodes;
					}
				}
			}
		}
	}
}
