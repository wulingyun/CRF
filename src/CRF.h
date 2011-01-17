#include <R.h>
#include <Rdefines.h>

/* Interfaces to R */

extern "C" {
	/* Decoding */
	SEXP Decode_Exact(SEXP _crf);
	SEXP Decode_Chain(SEXP _crf);
	SEXP Decode_Tree(SEXP _crf);
	SEXP Decode_LBP(SEXP _crf, SEXP _maxIter, SEXP _cutoff, SEXP _verbose);
	SEXP Decode_Sample(SEXP _crf, SEXP _samples);

	/* Inference */
	SEXP Infer_Exact(SEXP _crf);
	SEXP Infer_Chain(SEXP _crf);
	SEXP Infer_Tree(SEXP _crf);
	SEXP Infer_LBP(SEXP _crf, SEXP _maxIter, SEXP _cutoff, SEXP _verbose);
	SEXP Infer_Sample(SEXP _crf, SEXP _samples);

	/* Sampling */
	SEXP Sample_Exact(SEXP _crf, SEXP _size);
	SEXP Sample_Chain(SEXP _crf, SEXP _size);
	SEXP Sample_Tree(SEXP _crf, SEXP _size);
	SEXP Sample_Gibbs(SEXP _crf, SEXP _size, SEXP _burnIn, SEXP _start);

	/* Utils */
	SEXP Clamp_NodePot(SEXP _crfClamped);
}

/* CRF class */

class CRF {
public:
	SEXP _nNodes, _nEdges, _edges, _nStates, _maxState;
	int nNodes, nEdges, *edges, *nStates, maxState;

	SEXP _nAdj, _adjNodes, _adjEdges;
	int *nAdj, **adjNodes, **adjEdges;

	SEXP _nodePot, _edgePot;
	double *nodePot, *edgePot;

	CRF(SEXP _crf);
	~CRF();

	void Clamp_Reset(int *clamped, int *nodeMap, int nNodesNew, double *nodePotNew);

	/* BP functions */
	void TreeBP(double *messages_1, double *messages_2);
	void TreeBP_max(double *messages_1, double *messages_2);
	void LoopyBP(double *messages_1, double *messages_2, int maxIter, double cutoff, int verbose);
	void LoopyBP_max(double *messages_1, double *messages_2, int maxIter, double cutoff, int verbose);
	void Message2NodeBelief(double *messages_1, double *messages_2, double *nodeBel);
	void Message2EdgeBelief(double *messages_1, double *messages_2, double *nodeBel, double *edgeBel);
	void MaxOfMarginals(double *nodeBel, int *labels);
	void BetheFreeEnergy(double *nodeBel, double *edgeBel, double *logZ);

	/* Decoding methods */
	void Decode_Tree(int *labels);
	void Decode_LBP(int *labels, int maxIter, double cutoff, int verbose);
	/* Inference methods */
	void Infer_Tree(double *nodeBel, double *edgeBel, double *logZ);
	void Infer_LBP(double *nodeBel, double *edgeBel, double *logZ, int maxIter, double cutoff, int verbose);
	/* Sampling methods */
	void Sample_Tree(int size, int *samples);
	void Sample_LBP(int size, int *samples, int maxIter, double cutoff, int verbose);
};

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
