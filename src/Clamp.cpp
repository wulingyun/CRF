#include "CRF.h"

SEXP Clamp_Reset(SEXP _crf)
{
	CRFclamped crf(_crf);
	crf.Reset_NodePot();
	return(_crf);
}

CRFclamped::CRFclamped(SEXP _crf)
: CRF(_crf)
{
	_original = GetVar(_crf, "original");
	original.Set_Data(_original);

	PROTECT(_nodeId = AS_INTEGER(GetVar(_crf, "node.id")));
	PROTECT(_nodeMap = AS_INTEGER(GetVar(_crf, "node.map")));
	PROTECT(_edgeId = AS_INTEGER(GetVar(_crf, "edge.id")));
	PROTECT(_edgeMap = AS_INTEGER(GetVar(_crf, "edge.map")));
	nodeId = INTEGER_POINTER(_nodeId);
	nodeMap = INTEGER_POINTER(_nodeMap);
	edgeId = INTEGER_POINTER(_edgeId);
	edgeMap = INTEGER_POINTER(_edgeMap);

	PROTECT(_clamped = AS_INTEGER(GetVar(_crf, "clamped")));
	clamped = INTEGER_POINTER(_clamped);

	numProtect += 5;
}

void CRFclamped::Reset_NodePot()
{
	int e, n1, n2;

	for (int i = 0; i < original.nNodes; i++)
		if (nodeMap[i] > 0)
			for (int j = 0; j < original.nStates[i]; j++)
				NodePot(nodeMap[i]-1, j) = original.NodePot(i, j);

	for (int i = 0; i < original.nNodes; i++)
	{
		if (clamped[i])
		{
			for (int j = 0; j < original.nAdj[i]; j++)
			{
				e = original.AdjEdges(i, j);
				n1 = original.EdgesBegin(e);
				n2 = original.EdgesEnd(e);
				if (n1 == i && clamped[n2] == 0)
				{
					for (int k = 0; k < original.nStates[n2]; k++)
						NodePot(nodeMap[n2]-1, k) *= original.EdgePot(e, clamped[i]-1, k);
				}
				else if (n2 == i && clamped[n1] == 0)
				{
					for (int k = 0; k < original.nStates[n1]; k++)
						NodePot(nodeMap[n1]-1, k) *= original.EdgePot(e, k, clamped[i]-1);
				}
			}
		}
	}

	SetVar(_crf, "node.pot", _nodePot);
}
