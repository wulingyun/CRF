#include "CRF.h"

/* Residual BP */

void CRF::ResidualBP(int maxIter, double cutoff, int verbose, bool maximize)
{
  messages = (double ***) R_allocArray<double>(2, nEdges, maxState);
  double ***new_messages = (double ***) R_allocArray<double>(2, nEdges, maxState);
  double **priority = (double **) R_allocArray<double>(2, nEdges);

  for (int i = 0; i < nEdges; i++)
    for (int j = 0; j < maxState; j++)
    {
      messages[0][i][j] = new_messages[0][i][j] = 0;
      messages[1][i][j] = new_messages[1][i][j] = 0;
    }
    
  double *outgoing = (double *) R_alloc(maxState, sizeof(double));
  
  int s, r, e, n;

  for (int i = 0; i < nEdges; i++)
  {
    n = EdgesBegin(i);
    for (int j = 0; j < nStates[n]; j++)
      messages[0][i][j] = 1.0 / nStates[n];
    n = EdgesEnd(i);
    for (int j = 0; j < nStates[n]; j++)
      messages[1][i][j] = 1.0 / nStates[n];
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

      UpdateMessagePriority(s, r, e, messages, new_messages, priority);
    }
  }

  int d;
  double q, *msg, *new_msg;
  double difference = 0;
  for (int iter = 1; iter <= maxIter; iter++)
  {
    R_CheckUserInterrupt();

    q = 0;
    for (int i = 0; i < nEdges; i++)
    {
      Rprintf("   %d, %f, %f, %f\n", i, q, priority[0][i], priority[1][i]);
      if (priority[0][i] > q)
      {
        q = priority[0][i];
        e = i;
        d = 0;
      }
      if (priority[1][i] > q)
      {
        q = priority[1][i];
        e = i;
        d = 1;
      }
    }

    if (d == 0)
    {
      s = EdgesBegin(e);
      r = EdgesEnd(e);
      msg = messages[0][e];
      new_msg = messages[0][e];
      priority[d][e] = 0;
    }
    else
    {
      s = EdgesEnd(e);
      r = EdgesBegin(e);
      msg = messages[1][e];
      new_msg = messages[1][e];
      priority[d][e] = 0;
    }
    
    
    for (int i = 0; i < nStates[r]; i++)
    {
      msg[i] = new_msg[i];
    }

    GatherIncomingMessages(r, messages);

    Rprintf("Iter %d: %d, %d, %d\n", iter, s, r, e);
        
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
        
        UpdateMessagePriority(r, r0, e0, messages, new_messages, priority);
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

void CRF::UpdateMessagePriority(int s, int r, int e, double ***messages, double ***new_messages, double **priority)
{
  double *msg, *new_msg, *pri, res;
  if (EdgesBegin(e) == s)
  {
    msg = messages[0][e];
    new_msg = new_messages[0][e];
    pri = priority[0];
  }
  else
  {
    msg = messages[1][e];
    new_msg = new_messages[1][e];
    pri = priority[1];
  }
  
  pri[e] = 0;
  for (int i = 0; i < nStates[r]; i++)
  {
    res = fabs(msg[i] - new_msg[i]);
    if (res > pri[e]) pri[e] = res;
  }
}
