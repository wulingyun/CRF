#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include "misc.h"
using namespace std;

/* Interfaces to R */

extern "C" {
  /* DLL Init */
  void R_init_CRF(DllInfo *info);

	/* Decoding */
	SEXP Decode_Exact(SEXP _crf);
	SEXP Decode_Chain(SEXP _crf);
	SEXP Decode_Tree(SEXP _crf);
	SEXP Decode_Cutset(SEXP _crf, SEXP _engine, SEXP _start);
	SEXP Decode_Junction(SEXP _crf);
	SEXP Decode_Sample(SEXP _crf, SEXP _samples);
	SEXP Decode_LBP(SEXP _crf, SEXP _maxIter, SEXP _cutoff, SEXP _verbose);
	SEXP Decode_RBP(SEXP _crf, SEXP _maxIter, SEXP _cutoff, SEXP _verbose);
	SEXP Decode_TRBP(SEXP _crf, SEXP _maxIter, SEXP _cutoff, SEXP _verbose);
	SEXP Decode_Greedy(SEXP _crf, SEXP _restart, SEXP _start);
	SEXP Decode_ICM(SEXP _crf, SEXP _restart, SEXP _start);

	/* Inference */
	SEXP Infer_Exact(SEXP _crf);
	SEXP Infer_Chain(SEXP _crf);
	SEXP Infer_Tree(SEXP _crf);
	SEXP Infer_Cutset(SEXP _crf, SEXP _engine);
	SEXP Infer_Junction(SEXP _crf);
	SEXP Infer_Sample(SEXP _crf, SEXP _samples);
	SEXP Infer_LBP(SEXP _crf, SEXP _maxIter, SEXP _cutoff, SEXP _verbose, SEXP _maximize);
	SEXP Infer_RBP(SEXP _crf, SEXP _maxIter, SEXP _cutoff, SEXP _verbose, SEXP _maximize);
	SEXP Infer_TRBP(SEXP _crf, SEXP _maxIter, SEXP _cutoff, SEXP _verbose, SEXP _maximize);

	/* Sampling */
	SEXP Sample_Exact(SEXP _crf, SEXP _size);
	SEXP Sample_Chain(SEXP _crf, SEXP _size);
	SEXP Sample_Tree(SEXP _crf, SEXP _size);
	SEXP Sample_Cutset(SEXP _crf, SEXP _size, SEXP _engine);
	SEXP Sample_Junction(SEXP _crf, SEXP _size);
	SEXP Sample_Gibbs(SEXP _crf, SEXP _size, SEXP _burnIn, SEXP _start);

	/* Utils */
	SEXP Make_AdjInfo(SEXP _crf);
	SEXP MRF_Update(SEXP _crf);
	SEXP CRF_Update(SEXP _crf, SEXP _nodeFea, SEXP _edgeFea, SEXP _nodeExt, SEXP _edgeExt);
	SEXP MRF_Stat(SEXP _crf, SEXP _instances);
	SEXP MRF_NLL(SEXP _crf, SEXP _par, SEXP _instances, SEXP _infer, SEXP _env);
	SEXP CRF_NLL(SEXP _crf, SEXP _par, SEXP _instances, SEXP _nodeFea, SEXP _edgeFea, SEXP _nodeExt, SEXP _edgeExt, SEXP _infer, SEXP _env);
	SEXP Clamp_Reset(SEXP _crfClamped);
	SEXP Get_Potential(SEXP _crf, SEXP _configuration);
	SEXP Get_LogPotential(SEXP _crf, SEXP _configuration);
	SEXP Calc_Frequency(SEXP _v, SEXP _n);
}

class JunctionTree;

/* CRF class */

class CRF {
public:
	SEXP _crf;

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
	void Update_Pot(SEXP _nodeFea, SEXP _edgeFea, SEXP _nodeExt, SEXP _edgeExt);
	void Update_Pot();
	void Update_Pot_Finalize();
	double Get_Potential(int *configuration);
	double Get_LogPotential(int *configuration);
	void Normalize_NodePot();
	void Normalize_EdgePot();
	void Normalize_NodeBel();
	void Normalize_EdgeBel();
	void UB_Init();
	void UB_Clamp(int *clamped);
	double UB_Estimate();
	double UB_Estimate(int *clamped);

	/* BP functions */
	void TreeBP(bool maximize = false);
	void LoopyBP(int maxIter, double cutoff, int verbose, bool maximize = false);
	void ResidualBP(int maxIter, double cutoff, int verbose, bool maximize = false);

	void Messages2EdgeBel();
	void MaxOfMarginals();
	void BetheFreeEnergy();
	double *ComputeMessagesSum(int s, int r, int e, double *outgoing, double ***old_messages, double ***new_messages);
	double *ComputeMessagesMax(int s, int r, int e, double *outgoing, double ***old_messages, double ***new_messages);
	double UpdateMessagePriority(int s, int r, int e, double ***messages, double ***new_messages);
	void GatherIncomingMessages(int s, double ***old_messages);
	
	void TRBP(double *mu, double **scaleEdgePot, int maxIter, double cutoff, int verbose, bool maximize = false);
	void TRBP_Messages2EdgeBel(double *mu, double **scaleEdgePot);
	void TRBP_BetheFreeEnergy(double *mu);
	void TRBP_Init(double *mu, double **scaleEdgePot);

	/* Decoding methods */
	void Decode_Exact();
	void Decode_Chain();
	void Decode_Tree();
	void Decode_Junction();
	void Decode_Sample();
	void Decode_LBP(int maxIter, double cutoff, int verbose);
	void Decode_RBP(int maxIter, double cutoff, int verbose);
	void Decode_TRBP(int maxIter, double cutoff, int verbose);
	void Decode_Greedy(int restart = 0, int *start = 0);
	void Decode_ICM(int restart = 0, int *start = 0);

	/* Inference methods */
	void Infer_Exact();
	void Infer_Chain();
	void Infer_Tree();
	void Infer_Junction();
	void Infer_Sample();
	void Infer_LBP(int maxIter, double cutoff, int verbose, bool maximize = false);
	void Infer_RBP(int maxIter, double cutoff, int verbose, bool maximize = false);
	void Infer_TRBP(int maxIter, double cutoff, int verbose, bool maximize = false);

	/* Sampling methods */
	void Sample_Exact(int size = 0);
	void Sample_Chain(int size = 0);
	void Sample_Tree(int size = 0);
	void Sample_Junction(int size = 0);
	void Sample_Gibbs(int burnIn, int *start);
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

/* JunctionTree class */

class JunctionTree
{
public:
	CRF &original;
	int nNodes, nEdges, *nStates;

	int nClusters, **clusterNodes, *nClusterNodes, **clusterEdges, *nClusterEdges;
	int nSeperators, **seperatorNodes, *nSeperatorNodes, **seperators;
	int *nAdj, **adjClusters, **adjSeperators;
	int *nClusterStates, *nSeperatorStates;
	double **clusterBel, **seperatorBel;

	int cid, sid, *masks, *states;

	JunctionTree(CRF &crf);

	int States2Index(int nNodes, int *nodes, int *states);
	int *Index2States(int nNodes, int *nodes, int index, int *states);

	double &ClusterBel(int n, int *states);
	double &SeperatorBel(int n, int *states);

	void InitStateMasks(int c, int s = -1);
	void ResetClusterState();
	void ResetSeperatorState();
	bool NextClusterState();
	bool NextSeperatorState();

	void SendMessagesFromClusterSum(int c, int s);
	void SendMessagesFromClusterMax(int c, int s);
	void SendMessagesFromSeperator(int s, int c);
	void InitMessages();
	void Messages2NodeBel(bool maximize = false);
	void Messages2EdgeBel();

	void SendMessages(bool maximize = false);
	void Sample(int size);
};

inline double &JunctionTree::ClusterBel(int n, int *states)
{
	return clusterBel[n][States2Index(nClusterNodes[n], clusterNodes[n], states)];
}

inline double &JunctionTree::SeperatorBel(int n, int *states)
{
	return seperatorBel[n][States2Index(nSeperatorNodes[n], seperatorNodes[n], states)];
}
