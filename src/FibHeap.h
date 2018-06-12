
#ifndef __FIBHEAP_H
#define __FIBHEAP_H


#define OK      0
#define NOTOK   -1


class FibHeap;

class FibHeapNode {
  friend class FibHeap;
  FibHeapNode *m_left, *m_right, *m_parent, *m_child;
  short m_degree, m_mark, m_neg_infinity_flag;
  
protected:
  inline int  FHN_compare(FibHeapNode &RHS);
  inline void FHN_assign(FibHeapNode &RHS);
  
public:
  inline FibHeapNode();
  virtual ~FibHeapNode();
  
  virtual void operator =(FibHeapNode &RHS);
  virtual int  operator ==(FibHeapNode &RHS);
  virtual int  operator <(FibHeapNode &RHS);
  
  virtual void print();
  
};


class FibHeap {
  FibHeapNode *m_min_root;
  long m_num_nodes, m_num_trees, m_num_marked_nodes;
  
  int m_heap_ownership_flag;
  
public:
  FibHeap();
  virtual ~FibHeap();
  
  // The Standard Heap Operations
  
  void insert(FibHeapNode *new_node);
  void makeUnion(FibHeap *other_heap);
  
  inline FibHeapNode *minimum();
  FibHeapNode *extractMin();
  
  int decreaseKey(FibHeapNode *node, FibHeapNode &new_key);
  int remove(FibHeapNode *node);
  
  // Extra utility functions
  
  int  getHeapOwnership() { return m_heap_ownership_flag; };
  void setHeapOwnership() { m_heap_ownership_flag = 1; };
  void clearHeapOwnership() { m_heap_ownership_flag = 0; };
  
  long getNumNodes() { return m_num_nodes; };
  long getNumTrees() { return m_num_trees; };
  long getNumMarkedNodes() { return m_num_marked_nodes; };
  
  void print(FibHeapNode *tree = NULL, FibHeapNode *parent = NULL);
  
private:
  // Internal functions that help to implement the Standard Operations
  
  inline void exchange(FibHeapNode *&, FibHeapNode *&);
  void consolidate();
  void link(FibHeapNode *, FibHeapNode *);
  void addToRootList(FibHeapNode *);
  void cut(FibHeapNode *, FibHeapNode *);
  void cascadingCut(FibHeapNode *);
};


FibHeapNode::FibHeapNode()
{
  m_left = m_right = m_parent = m_child = NULL;
  m_degree = m_mark = m_neg_infinity_flag = 0;
}

//=========================================================
// FHN_assign()
//
// To be used as first step of an assignment operator in a
// derived class.  The derived class will handle assignment
// of key value, and this function handles copy of the
// m_neg_infinity_flag (which overrides the key value if it is
// set).
//=========================================================

void FibHeapNode::FHN_assign(FibHeapNode &RHS)
{
  m_neg_infinity_flag = RHS.m_neg_infinity_flag;
}

//=========================================================
// FHN_compare()
//
// To be used as the first step of ALL comparators in a
// derived class.
//
// Compares the relative state of the two neg. infinity
// flags.  Note that 'this' is the m_left hand side.  If
// LHS neg. infinity is set, then it will be less than (-1)
// the RHS unless RHS neg. infinity flag is also set.
// Only if function returns 0 should the key comparison
// defined in the derived class be performed, e.g.
//
// For ==, if zero returned, then compare keys
//	   if non-zero X returned, then return 0
// For <,  if zero returned, then compare keys
//         if non-zero X returned, then return X<0?1:0
// For >,  if zero returned, then compare keys
//         if non-zero X returned, then return X>0?1:0    
//=========================================================

int  FibHeapNode::FHN_compare(FibHeapNode &RHS)
{
  if (m_neg_infinity_flag) return RHS.m_neg_infinity_flag ? 0 : -1;
  return RHS.m_neg_infinity_flag ? 1 : 0; 
}


#endif // __FIBHEAP_H
