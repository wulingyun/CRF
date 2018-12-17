#include "CRF.h"

SEXP Decode_Cutset(SEXP _crf, SEXP _engine, SEXP _start)
{
	CRFclamped crf(_crf);
	crf.Init_Labels();
	crf.Init_NodeBel();
	crf.original.Init_Labels();
	
	PROTECT(_start = AS_INTEGER(_start));
	int *start = INTEGER_POINTER(_start);

	crf.Decode_Cutset(INTEGER_POINTER(AS_INTEGER(_engine))[0], start);

	UNPROTECT(1);

	return(crf.original._labels);
}

void CRFclamped::Decode_Cutset(int engine, int *start)
{
	original.UB_Init();
	original.UB_Clamp(clamped);

	int *y = (int *) R_alloc(original.nNodes, sizeof(int));
	double max;
	if (start)
		for (int i = 0; i < original.nNodes; i++)
			y[i] = start[i] - 1;
	else
		for (int i = 0; i < original.nNodes; i++)
		{
			max = -1;
			for (int j = 0; j < original.nStates[i]; j++)
				if (max < original.NodePot(i,j))
				{
					max = original.NodePot(i,j);
					y[i] = j;
				}
		}
	double maxPot = original.Get_Potential(y);
	for (int i = 0; i < original.nNodes; i++)
		original.labels[i] = y[i] + 1;

	for (int i = 0; i < original.nNodes; i++)
	{
		if (clamped[i] > 0)
		{
			clamped[i] = 1;
			y[i] = 0;
		}
		else
		{
			clamped[i] = 0;
			y[i] = -1;
		}
	}

	double pot;
	int index;
	while (1)
	{
		R_CheckUserInterrupt();

		if (original.UB_Estimate(clamped) > maxPot)
		{
			/* Reset node potentials */
			Reset_NodePot();

			/* Decode clamped CRF */
			switch (engine)
			{
			case 0:
				break;
			case 1:
				Decode_Exact();
				break;
			case 2:
				Decode_Chain();
				break;
			case 3:
				Decode_Tree();
				break;
			default:
				Decode_Tree();
				break;
			}

			/* Map results back */
			for (int i = 0; i < nNodes; i++)
				y[nodeId[i]-1] = labels[i] - 1;

			/* Calculate potential */
			pot = original.Get_Potential(y);

			/* Record the best potentials */
			if (pot > maxPot)
			{
				maxPot = pot;
				for (int i = 0; i < original.nNodes; i++)
					original.labels[i] = y[i] + 1;
			}
		}

		/* Next configuration */
		for (index = 0; index < original.nNodes; index++)
		{
			if (clamped[index] == 0)
				continue;
			clamped[index]++;
			y[index]++;
			if (y[index] < original.nStates[index])
				break;
			else
			{
				clamped[index] = 1;
				y[index] = 0;
			}
		}

		if (index == original.nNodes)
			break;
	}
}
