#include "CRF.h"

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

	UNPROTECT(2);
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

	UNPROTECT(2);
	return(_potential);
}
