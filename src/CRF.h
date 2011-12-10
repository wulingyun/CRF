#include <R.h>
#include <Rdefines.h>

/* Interfaces to R */

extern "C" {
	/* Decoding */
	SEXP Decode_Exact(SEXP _crf);
	SEXP Decode_Chain(SEXP _crf);
	SEXP Decode_Tree(SEXP _crf);
	SEXP Decode_Cutset(SEXP _crf, SEXP _engine, SEXP _start);
	SEXP Decode_Sample(SEXP _crf, SEXP _samples);
	SEXP Decode_LBP(SEXP _crf, SEXP _maxIter, SEXP _cutoff, SEXP _verbose);
	SEXP Decode_TRBP(SEXP _crf, SEXP _maxIter, SEXP _cutoff, SEXP _verbose);
	SEXP Decode_Greedy(SEXP _crf, SEXP _restart, SEXP _start);
	SEXP Decode_ICM(SEXP _crf, SEXP _restart, SEXP _start);

	/* Inference */
	SEXP Infer_Exact(SEXP _crf);
	SEXP Infer_Chain(SEXP _crf);
	SEXP Infer_Tree(SEXP _crf);
	SEXP Infer_Cutset(SEXP _crf, SEXP _engine);
	SEXP Infer_Sample(SEXP _crf, SEXP _samples);
	SEXP Infer_LBP(SEXP _crf, SEXP _maxIter, SEXP _cutoff, SEXP _verbose);
	SEXP Infer_TRBP(SEXP _crf, SEXP _maxIter, SEXP _cutoff, SEXP _verbose);

	/* Sampling */
	SEXP Sample_Exact(SEXP _crf, SEXP _size);
	SEXP Sample_Chain(SEXP _crf, SEXP _size);
	SEXP Sample_Tree(SEXP _crf, SEXP _size);
	SEXP Sample_Cutset(SEXP _crf, SEXP _size, SEXP _engine);
	SEXP Sample_Gibbs(SEXP _crf, SEXP _size, SEXP _burnIn, SEXP _start);

	/* Utils */
	SEXP Clamp_NodePot(SEXP _crfClamped);
	SEXP Get_Potential(SEXP _crf, SEXP _configuration);
	SEXP Get_LogPotential(SEXP _crf, SEXP _configuration);
}

/* CRF class */

class CRF {
public:
	SEXP _nNodes, _nEdges, _edges, _nStates, _maxState;
	int nNodes, nEdges, *edges, *nStates, maxState;

	SEXP _nAdj, _adjNodes, _adjEdges;
	int *nAdj, **adjNodes, **adjEdges;

	SEXP _nodePot, _edgePot;
	double *nodePot, **edgePot;
	int *nEdgeStates;

	SEXP _labels;
	int *labels;

	SEXP _nodeBel, _edgeBel, _logZ, _belief;
	double *nodeBel, **edgeBel, *logZ;

	SEXP _samples;
	int *samples, nSamples;

	double *maxNodePot, *maxEdgePot, unclampedUB;
	double ***messages;

	int numProtect;

	CRF();
	CRF(SEXP _crf);
	~CRF();

	/* Initialize results */
	void Init_Labels();
	void Init_NodeBel();
	void Init_EdgeBel();
	void Init_LogZ();
	void Init_Belief();
	void Init_Samples(int size);
	void Init_Samples(SEXP _size);

	/* Set members */
	void Set_Data(SEXP _crf);
	void Set_Samples(SEXP _otherSamples);

	int AdjNodes(int n, int i);
	int AdjEdges(int n, int i);
	int EdgesBegin(int e);
	int EdgesEnd(int e);
	double &NodePot(int n, int s);
	double &EdgePot(int e, int s1, int s2);
	double &NodeBel(int n, int s);
	double &EdgeBel(int e, int s1, int s2);
	int &Samples(int i, int n);

	/* Utils */
	double Get_Potential(int *configuration);
	double Get_LogPotential(int *configuration);
	void UB_Init();
	void UB_Clamp(int *clamped);
	double UB_Estimate();
	double UB_Estimate(int *clamped);

	/* BP functions */
	void TreeBP(bool maximize = false);
	void LoopyBP(int maxIter, double cutoff, int verbose, bool maximize = false);
	void MessagesInit();
	void Messages2NodeBel();
	void Messages2EdgeBel();
	void MaxOfMarginals();
	void BetheFreeEnergy();
	void TRBP(double *mu, double **scaleEdgePot, int maxIter, double cutoff, int verbose, bool maximize = false);
	void TRBP_Messages2NodeBel(double *mu);
	void TRBP_Messages2EdgeBel(double *mu, double **scaleEdgePot);
	void TRBP_BetheFreeEnergy(double *mu);
	void TRBP_Weights(double *mu);
	void TRBP_ScaleEdgePot(double *mu, double **scaleEdgePot);
	void TRBP_MinSpanTree(int *tree, double *costs);

	/* Decoding methods */
	void Decode_Exact();
	void Decode_Chain();
	void Decode_Tree();
	void Decode_Sample();
	void Decode_LBP(int maxIter, double cutoff, int verbose);
	void Decode_TRBP(int maxIter, double cutoff, int verbose);
	void Decode_Greedy(int restart = 0, int *start = 0);
	void Decode_ICM(int restart = 0, int *start = 0);

	/* Inference methods */
	void Infer_Exact();
	void Infer_Chain();
	void Infer_Tree();
	void Infer_Sample();
	void Infer_LBP(int maxIter, double cutoff, int verbose);
	void Infer_TRBP(int maxIter, double cutoff, int verbose);

	/* Sampling methods */
	void Sample_Exact(int size = 0);
	void Sample_Chain(int size = 0);
	void Sample_Tree(int size = 0);
	void Sample_Gibbs(int burnIn, int *start, int size = 0);
};

class CRFclamped: public CRF {
public:
	SEXP _original;
	CRF original;

	SEXP _clamped, _nodeId, _nodeMap, _edgeId, _edgeMap;
	int *clamped, *nodeId, *nodeMap, *edgeId, *edgeMap;

	CRFclamped(SEXP _crf);

	void Reset_NodePot();

	/* Decoding methods */
	void Decode_Cutset(int engine = -1, int *start = 0);

	/* Inference methods */
	void Infer_Cutset(int engine = -1);

	/* Sampling methods */
	void Sample_Cutset(int size = 0, int engine = -1);
};

/* inline functions */

inline int CRF::AdjNodes(int n, int i)
{
	return adjNodes[n][i] - 1;
}

inline int CRF::AdjEdges(int n, int i)
{
	return adjEdges[n][i] - 1;
}

inline int CRF::EdgesBegin(int e)
{
	return edges[e] - 1;
}

inline int CRF::EdgesEnd(int e)
{
	return edges[e + nEdges] - 1;
}

inline double &CRF::NodePot(int n, int s)
{
	return nodePot[n + nNodes * s];
}

inline double &CRF::EdgePot(int e, int s1, int s2)
{
	return edgePot[e][s1 + nStates[EdgesBegin(e)] * s2];
}

inline double &CRF::NodeBel(int n, int s)
{
	return nodeBel[n + nNodes * s];
}

inline double &CRF::EdgeBel(int e, int s1, int s2)
{
	return edgeBel[e][s1 + nStates[EdgesBegin(e)] * s2];
}

inline int &CRF::Samples(int i, int n)
{
	return samples[i + nSamples * n];
}


/* initialize the list */
template <class T>
void setValues(SEXP r, T *c, T v)
{
	for (int i = 0; i < length(r); i++)
		c[i] = v;
};

/* get/set the list element */
SEXP getListElement(SEXP list, const char *tag);
void setListElement(SEXP list, int i, const char *tag, SEXP value);

/* set dim of array */
void setDim2(SEXP array, int x1, int x2);
void setDim3(SEXP array, int x1, int x2, int x3);

/* sample from discret distribution */
int sample(int n, double *prob);

/* swap variables */
template <class T>
void swap(T &a, T &b)
{
	T temp = a;
	a = b;
	b = temp;
};

/* allocate multidimensional array */
template <class T, int n>
void **allocArray(int dim[n])
{
	int size1, size2;
	void **array, **sub1, **sub2, **tmp;
	size1 = dim[0];
	array = sub1 = (void **) R_alloc(size1, sizeof(void *));
	for (int i = 1; i < n-1; i++)
	{
		size2 = size1 * dim[i];
		tmp = sub2 = (void **) R_alloc(size2, sizeof(void *));
		for (int j = 0; j < size1; j++)
		{
			sub1[j] = (void *) tmp;
			tmp += dim[i];
		}
		size1 = size2;
		sub1 = sub2;
	}
	T *block = (T *) R_alloc(size1 * dim[n-1], sizeof(T));
	for (int i = 0; i < size1; i++)
	{
		sub1[i] = block;
		block += dim[n-1];
	}
	return array;
};

template <class T>
T **allocArray2(int dim1, int *dim2)
{
	T *block, **array;
	int array_size;
	array_size = 0;
	for (int i = 0; i < dim1; i++)
		array_size += dim2[i];
	block = (T *) R_alloc(array_size, sizeof(T));
	array = (T **) R_alloc(dim1, sizeof(T *));
	for (int i = 0; i < dim1; i++)
	{
		array[i] = block;
		block += dim2[i];
	}
	return array;
};
