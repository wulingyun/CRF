#include "CRF.h"
#include "FibHeap.h"
#include <Rmath.h>

// for using FibHeap

class HeapNode: public FibHeapNode
{
  double m_priority;
  int m_dir, m_index;
  
public:
  HeapNode() : FibHeapNode() { m_priority = 0; };
  
  virtual void operator =(HeapNode& RHS);
  virtual int  operator ==(HeapNode& RHS);
  virtual int  operator <(HeapNode& RHS);
  
  virtual void operator =(double key);
  
  double getKeyValue() { return m_priority; }
  void setKeyValue(double key) { m_priority = key; }
  
  int getDirValue() { return m_dir; }
  void setDirValue(int i) { m_dir = i; }
  int getIndexValue() { return m_index; }
  void setIndexValue(int i) { m_index = i; }
};

void HeapNode::operator =(double priority)
{
  HeapNode temp;
  temp.m_priority = m_priority = priority;
  FHN_assign(temp);
}

void HeapNode::operator =(HeapNode& RHS)
{
  FHN_assign(RHS);
  m_priority = RHS.m_priority;
}

int  HeapNode::operator ==(HeapNode& RHS)
{
  if (FHN_compare(RHS)) return 0;
  return m_priority == RHS.m_priority ? 1 : 0;
}

int  HeapNode::operator <(HeapNode& RHS)
{
  int X;
  
  if ((X=FHN_compare(RHS)) != 0) return X < 0 ? 1 : 0;
  return m_priority < RHS.m_priority ? 1 : 0;
}

/* Residual BP */

void CRF::ResidualBP(int maxIter, double cutoff, int verbose, bool maximize)
{
  messages = (double ***) R_allocArray<double>(2, nEdges, maxState);
  double ***new_messages = (double ***) R_allocArray<double>(2, nEdges, maxState);
  HeapNode *priority[2];
  FibHeap priority_heap;
  int d;
  
  priority[0] = new HeapNode [nEdges];
  priority[1] = new HeapNode [nEdges];
  
  for (int i = 0; i < nEdges; i++)
    for (int j = 0; j < maxState; j++)
    {
      messages[0][i][j] = new_messages[0][i][j] = 0;
      messages[1][i][j] = new_messages[1][i][j] = 0;
    }

  double *outgoing = (double *) R_alloc(maxState, sizeof(double));
  
  int s, r, e, n;

  double sumMsg;
  for (int i = 0; i < nEdges; i++)
  {
    sumMsg = 0;
    n = EdgesBegin(i);
    for (int j = 0; j < nStates[n]; j++)
      sumMsg += messages[0][i][j] = runif(0, 1);
    for (int j = 0; j < nStates[n]; j++)
      messages[0][i][j] /= sumMsg;

    sumMsg = 0;
    n = EdgesEnd(i);
    for (int j = 0; j < nStates[n]; j++)
      sumMsg += messages[1][i][j] = runif(0, 1);
    for (int j = 0; j < nStates[n]; j++)
      messages[1][i][j] /= sumMsg;
  }

  /* compute possible message update */
  
  for (s = 0; s < nNodes; s++)
  {
    GatherIncomingMessages(s, messages);
    
    for (int i = 0; i < nAdj[s]; i++)
    {
      r = AdjNodes(s, i);
      e = AdjEdges(s, i);
      
      if (maximize)
        ComputeMessagesMax(s, r, e, outgoing, messages, new_messages);
      else
        ComputeMessagesSum(s, r, e, outgoing, messages, new_messages);

      d = EdgesBegin(e) == s ? 1 : 0;
      priority[d][e].setKeyValue(UpdateMessagePriority(s, r, e, messages, new_messages));
      priority[d][e].setDirValue(d);
      priority[d][e].setIndexValue(e);
      priority_heap.insert(&priority[d][e]);
    }
  }

  double *msg, *new_msg;
  double difference = 0;
  for (int iter = 1; iter <= maxIter; iter++)
  {
    R_CheckUserInterrupt();
    
    for (int iterI = 0; iterI < nEdges; iterI++)
    {
      HeapNode *p = (HeapNode *) priority_heap.extractMin();
      p->setKeyValue(0);
      priority_heap.insert(p);
      
      d = p->getDirValue();
      e = p->getIndexValue();
      
      if (d == 0)
      {
        s = EdgesEnd(e);
        r = EdgesBegin(e);
        msg = messages[0][e];
        new_msg = new_messages[0][e];
      }
      else
      {
        s = EdgesBegin(e);
        r = EdgesEnd(e);
        msg = messages[1][e];
        new_msg = new_messages[1][e];
      }
      
      
      for (int i = 0; i < nStates[r]; i++)
      {
        msg[i] = new_msg[i];
      }
      
      GatherIncomingMessages(r, messages);
      
      int r0, e0;    
      for (int i = 0; i < nAdj[r]; i++)
      {
        r0 = AdjNodes(r, i);
        e0 = AdjEdges(r, i);
        if (r0 != s)
        {
          if (maximize)
            ComputeMessagesMax(r, r0, e0, outgoing, messages, new_messages);
          else
            ComputeMessagesSum(r, r0, e0, outgoing, messages, new_messages);
          
          d = EdgesBegin(e0) == r ? 1 : 0;
          priority_heap.remove(&priority[d][e0]);
          priority[d][e0].setKeyValue(UpdateMessagePriority(r, r0, e0, messages, new_messages));
          priority_heap.insert(&priority[d][e0]);
        }
      }
    }

    difference = 0;
    for (int i = 0; i < nEdges; i++)
    {
      for (int j = 0; j < maxState; j++)
      {
        difference += fabs(messages[0][i][j] - new_messages[0][i][j]);
        difference += fabs(messages[1][i][j] - new_messages[1][i][j]);
      }
    }

    if (verbose)
      Rprintf("ResidualBP: Iteration %d, Difference = %f\n", iter, difference);
    if (difference <= cutoff)
      break;
  }
  
  if (difference > cutoff)
    warning("Residual BP did not converge in %d iterations! (diff = %f)", maxIter, difference);
}

double CRF::UpdateMessagePriority(int s, int r, int e, double ***messages, double ***new_messages)
{
  double *msg, *new_msg, priority, res;
  if (EdgesBegin(e) == r)
  {
    msg = messages[0][e];
    new_msg = new_messages[0][e];
  }
  else
  {
    msg = messages[1][e];
    new_msg = new_messages[1][e];
  }
  
  priority = 0;
  for (int i = 0; i < nStates[r]; i++)
  {
    res = fabs(msg[i] - new_msg[i]);
    if (res > priority) priority = res;
  }
  
  return(-priority);
}
