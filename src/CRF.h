#include <R.h>
#include <Rdefines.h>


typedef struct __CRFinfo {
	SEXP _nNodes, _nEdges, _edges, _nStates, _maxState;
	int nNodes, nEdges, *edges, *nStates, maxState;

	SEXP _nAdj, _adjNodes, _adjEdges;
	int *nAdj, **adjNodes, **adjEdges;

	SEXP _nodePot, _edgePot;
	double *nodePot, *edgePot;
} CRFinfo;

void openCRF(CRFinfo *crf, SEXP _crf);
void closeCRF(CRFinfo *crf);

void TreeBP(CRFinfo *crf, double *messages_1, double *messages_2);
void TreeBP_max(CRFinfo *crf, double *messages_1, double *messages_2);
void LoopyBP(CRFinfo *crf, double *messages_1, double *messages_2, int maxIter, double cutoff, int verbose);
void LoopyBP_max(CRFinfo *crf, double *messages_1, double *messages_2, int maxIter, double cutoff, int verbose);
void Message2NodeBelief(CRFinfo *crf, double *messages_1, double *messages_2, double *nodeBel);
void Message2EdgeBelief(CRFinfo *crf, double *messages_1, double *messages_2, double *nodeBel, double *edgeBel);

void _Decode_Tree(CRFinfo *crf, int *labels);
void _Decode_LBP(CRFinfo *crf, int *labels, int maxIter, double cutoff, int verbose);
void _Infer_Tree(CRFinfo *crf, double *nodeBel, double *edgeBel, double *logZ);
void _Infer_LBP(CRFinfo *crf, double *nodeBel, double *edgeBel, double *logZ, int maxIter, double cutoff, int verbose);
void _Sample_Tree(CRFinfo *ctf, int size, int *samples);
void _Sample_LBP(CRFinfo *ctf, int size, int *samples, int maxIter, double cutoff, int verbose);

/* initialize the list */
#define setValues(r, c, v) for (int i = 0; i < length(r); i++) c[i] = v

/* get/set the list element */
SEXP getListElement(SEXP list, const char *tag);
void setListElement(SEXP list, int i, const char *tag, SEXP value);

/* set dim of array */
void setDim2(SEXP array, int x1, int x2);
void setDim3(SEXP array, int x1, int x2, int x3);

/* sample from discret distribution */
int sample(int n, double *prob);
