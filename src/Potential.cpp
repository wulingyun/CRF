#include "CRF.h"

SEXP Get_Potential(SEXP _crf, SEXP _configuration)
{
	CRF crf(_crf);

	PROTECT(_configuration = AS_INTEGER(_configuration));
	int *configuration = INTEGER_POINTER(_configuration);

	SEXP _potential;
	PROTECT(_potential = NEW_NUMERIC(1));
	double *potential = NUMERIC_POINTER(_potential);

	int *y = (int *) R_alloc(crf.nNodes, sizeof(int));
	for (int i = 0; i < crf.nNodes; i++)
		y[i] = configuration[i] - 1;

	*potential = crf.Get_Potential(y);

	UNPROTECT_PTR(_configuration);
	UNPROTECT_PTR(_potential);
	return(_potential);
}

SEXP Get_LogPotential(SEXP _crf, SEXP _configuration)
{
	CRF crf(_crf);

	PROTECT(_configuration = AS_INTEGER(_configuration));
	int *configuration = INTEGER_POINTER(_configuration);

	SEXP _potential;
	PROTECT(_potential = NEW_NUMERIC(1));
	double *potential = NUMERIC_POINTER(_potential);

	int *y = (int *) R_alloc(crf.nNodes, sizeof(int));
	for (int i = 0; i < crf.nNodes; i++)
		y[i] = configuration[i] - 1;

	*potential = crf.Get_LogPotential(y);

	UNPROTECT_PTR(_configuration);
	UNPROTECT_PTR(_potential);
	return(_potential);
}

double CRF::Get_Potential(int *configuration)
{
	double potential = 1;

	/* Node potentials */
	for (int i = 0; i < nNodes; i++)
		potential *= nodePot[i + nNodes * configuration[i]];

	/* Edge potentials */
	for (int i = 0; i < nEdges; i++)
		potential *= edgePot[configuration[edges[i]-1] + maxState * (configuration[edges[i+nEdges]-1] + maxState * i)];

	return(potential);
}

double CRF::Get_LogPotential(int *configuration)
{
	double potential = 0;

	/* Node potentials */
	for (int i = 0; i < nNodes; i++)
		potential += log(nodePot[i + nNodes * configuration[i]]);

	/* Edge potentials */
	for (int i = 0; i < nEdges; i++)
		potential += log(edgePot[configuration[edges[i]-1] + maxState * (configuration[edges[i+nEdges]-1] + maxState * i)]);

	return(potential);
}

void CRF::UB_Init()
{
	maxNodePot = (double *) R_alloc(nNodes, sizeof(double));
	maxEdgePot = (double *) R_alloc(nEdges, sizeof(double));

	/* Node potentials */
	for (int i = 0; i < nNodes; i++)
	{
		maxNodePot[i] = 0;
		for (int j = 0; j < nStates[i]; j++)
			if (maxNodePot[i] < NodePot(i, j))
				maxNodePot[i] = NodePot(i, j);
	}

	/* Edge potentials */
	int n1, n2;
	for (int i = 0; i < nEdges; i++)
	{
		maxEdgePot[i] = 0;
		n1 = Edges(i,0) - 1;
		n2 = Edges(i,1) - 1;
		for (int j = 0; j < nStates[n1]; j++)
			for (int k = 0; k < nStates[n2]; k++)
				if (maxEdgePot[i] < EdgePot(j, k, i))
					maxEdgePot[i] = EdgePot(j, k, i);
	}
}

void CRF::UB_Clamp(int *clamped)
{
	unclampedUB = 1;
	int s1, s2;

	/* Node potentials */
	for (int i = 0; i < nNodes; i++)
		if (clamped[i] <= 0)
			unclampedUB *= maxNodePot[i];

	/* Edge potentials */
	for (int i = 0; i < nEdges; i++)
	{
		s1 = clamped[Edges(i,0) - 1] - 1;
		s2 = clamped[Edges(i,1) - 1] - 1;
		if (s1 < 0 || s2 < 0)
			unclampedUB *= maxEdgePot[i];
	}
}

double CRF::UB_Estimate()
{
	double potential = 1;

	/* Node potentials */
	for (int i = 0; i < nNodes; i++)
		potential *= maxNodePot[i];

	/* Edge potentials */
	for (int i = 0; i < nEdges; i++)
		potential *= maxEdgePot[i];

	return(potential);
}

double CRF::UB_Estimate(int *clamped)
{
	double potential = unclampedUB;
	int s1, s2;

	/* Node potentials */
	for (int i = 0; i < nNodes; i++)
		if (clamped[i] > 0)
			potential *= NodePot(i, clamped[i]-1);

	/* Edge potentials */
	for (int i = 0; i < nEdges; i++)
	{
		s1 = clamped[Edges(i,0) - 1] - 1;
		s2 = clamped[Edges(i,1) - 1] - 1;
		if (s1 >= 0 && s2 >= 0)
			potential *= EdgePot(s1, s2, i);
	}

	return(potential);
}
