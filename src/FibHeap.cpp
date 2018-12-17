//***************************************************************************
// This Fibonacci heap implementation is Copym_right (c) 1996 by John Boyer.
// See the header file for free usage information.
//***************************************************************************

//***************************************************************************
// The classes in this package are designed to allow the package user
// to quickly and easily develop applications that require a heap data
// structure.  Using amortized analysis, the asymptotically fastest heap
// data structure is the Fibonacci heap.  The constants are a little
// high so the real speed gain will not be seen until larger data sets
// are required, but in most cases, if the data set is small, then the
// run-time will be neglible anyway.
//
// To use this heap class you need do only two things.  First, subclass
// the FibHeapNode class to create the class of objects that you'd
// like to store in a heap.  Second, create an instance of the FibHeap
// class, which can then be used to Insert(), extractMin(), etc.,
// instances of your FibHeapNode subclass.  Notice that you don't need
// to create a subclass of the FibHeap class.
//
// The application-specific data object that you'd like to store in a heap
// will have a key value.  In the class that you derive from FibHeapNode,
// you will need to define the key structure then provide assignment (=),
// equality (==) and less-than operators and a destructor.  These functions
// are declared virtual so that the code in the FibHeap class can compare,
// assign and destroy your objects by calling on your code.
//
// The overloaded operators in your defined class MUST call functions in
// the Fibonacci heap node class first.  For assignment, the function
// FHN_Assign() should be called before code that deals with the copy of
// the key value.  For comparison operators, the function FHN_Cmp() should
// appear first.  If it returns 0, then keys can be compared as normal.
// The following indicates what the three most common operators must do
// based on the return value of FHN_Cmp() 
//
// For ==, if zero returned, then compare keys
//	   if non-zero X returned, then return 0
// For <,  if zero returned, then compare keys
//         if non-zero X returned, then return X<0?1:0
// For >,  if zero returned, then compare keys
//         if non-zero X returned, then return X>0?1:0   
//***************************************************************************


#include <stdlib.h>
//#include <iostream.h>
#include <stdio.h>
//#include <conio.h>
#include <R.h>

#include "FibHeap.h"


//=========================================================
// FibHeapNode Destructor
//
// Body is empty, but declaration is required in order to
// force virtual.  This will ensure that FibHeap class
// calls derived class destructors.
//=========================================================

FibHeapNode::~FibHeapNode()
{
}

//========================================================================
// We do, on occasion, compare and assign objects of type FibHeapNode, but
// only when the m_neg_infinity_flag is set.  See for example FibHeap::Delete().
//
// Also, these functions exemplify what a derived class should do.
//========================================================================

void FibHeapNode::operator =(FibHeapNode &RHS)
{
  FHN_assign(RHS);
  // Key assignment goes here in derived classes
}

int  FibHeapNode::operator ==(FibHeapNode &RHS)
{
  if (FHN_compare(RHS)) return 0;
  // Key compare goes here in derived classes
  return 1;
}

int  FibHeapNode::operator <(FibHeapNode &RHS)
{
  int i;
  
  if ((i=FHN_compare(RHS)) != 0) return i < 0 ? 1 : 0;
  // Key compare goes here in derived classes
  return 0;
}

//=========================================================
// Print()
//=========================================================

void FibHeapNode::print()
{
  //	if (m_neg_infinity_flag) cout << "-inf.";
}

//***************************************************************************
//===========================================================================
// FibHeap Constructor
//===========================================================================
//***************************************************************************

FibHeap::FibHeap()
{
  m_min_root = NULL;
  m_num_nodes = m_num_trees = m_num_marked_nodes = 0;
  clearHeapOwnership();
}

//===========================================================================
// FibHeap Destructor
//===========================================================================

FibHeap::~FibHeap()
{
  FibHeapNode *temp;
  
  if (getHeapOwnership()) {
    while (m_min_root != NULL) {
      temp = extractMin();
      delete temp;
    }
  }
}

//===========================================================================
// Insert() - O(1) actual; O(2) amortized
//
// I am using O(2) here to indicate that although Insert() is
// constant time, its amortized rating is more costly because some
// of the work of inserting is done by other operations such as
// extractMin(), which is where tree-balancing occurs.
//
// The m_child pointer is deliberately not set to NULL because Insert()
// is also used internally to help put whole trees onto the root list.
//===========================================================================

void FibHeap::insert(FibHeapNode *NewNode)
{
  if (NewNode == NULL) return;
  
  // If the heap is currently empty, then new node becomes singleton
  // circular root list
  
  if (m_min_root == NULL) {
    m_min_root = NewNode->m_left = NewNode->m_right = NewNode;
  }
  else {
    // Pointers from NewNode set to insert between m_min_root and m_min_root->m_right
    
    NewNode->m_right = m_min_root->m_right;
    NewNode->m_left = m_min_root;
    
    // Set Pointers to NewNode  
    
    NewNode->m_left->m_right = NewNode;
    NewNode->m_right->m_left = NewNode;
    
    // The new node becomes new m_min_root if it is less than current m_min_root
    
    if (*NewNode < *m_min_root) m_min_root = NewNode;
  }
  
  // We have one more node in the heap, and it is a tree on the root list
  
  m_num_nodes++;
  
  m_num_trees++;
  NewNode->m_parent = NULL;
}

//===========================================================================
// Union() - O(1) actual; O(1) amortized
//===========================================================================

void FibHeap::makeUnion(FibHeap *OtherHeap)
{
  FibHeapNode *Min1, *Min2, *Next1, *Next2;
  
  if (OtherHeap == NULL || OtherHeap->m_min_root == NULL) return;
  
  // We join the two circular lists by cutting each list between its
  // min node and the node after the min.  This code just pulls those
  // nodes into temporary variables so we don't get lost as changes
  // are made.
  
  Min1 = m_min_root;
  Min2 = OtherHeap->m_min_root;
  Next1 = Min1->m_right;
  Next2 = Min2->m_right;
  
  // To join the two circles, we join the minimum nodes to the next
  // nodes on the opposite chains.  Conceptually, it looks like the way
  // two bubbles join to form one larger bubble.  They meet at one point
  // of contact, then expand out to make the bigger circle.
  
  Min1->m_right = Next2;
  Next2->m_left = Min1;
  Min2->m_right = Next1;
  Next1->m_left = Min2;
  
  // Choose the new minimum for the heap
  
  if (*Min2 < *Min1) m_min_root = Min2;
  
  // Set the amortized analysis statistics and size of the new heap
  
  m_num_nodes += OtherHeap->m_num_nodes;
  m_num_marked_nodes += OtherHeap->m_num_marked_nodes;
  m_num_trees += OtherHeap->m_num_trees;
  
  // Complete the union by setting the other heap to emptiness
  // then destroying it
  
  OtherHeap->m_min_root  = NULL;
  OtherHeap->m_num_nodes =
    OtherHeap->m_num_trees =
    OtherHeap->m_num_marked_nodes = 0;
  
  delete OtherHeap;
}

//===========================================================================
// Minimum - O(1) actual; O(1) amortized
//===========================================================================

FibHeapNode *FibHeap::minimum()
{
  return m_min_root;
}

//===========================================================================
// extractMin() - O(n) worst-case actual; O(lg n) amortized
//===========================================================================

FibHeapNode *FibHeap::extractMin()
{
  FibHeapNode *Result;
  FibHeap *m_childHeap = NULL;
  
  // Remove minimum node and set m_min_root to next node
  
  if ((Result = minimum()) == NULL) return NULL;
  
  m_min_root = Result->m_right;
  Result->m_right->m_left = Result->m_left;
  Result->m_left->m_right = Result->m_right;
  Result->m_left = Result->m_right = NULL;
  
  m_num_nodes --;
  if (Result->m_mark) {
    m_num_marked_nodes --;
    Result->m_mark = 0;
  }
  Result->m_degree = 0;
  
  // Attach m_child list of Minimum node to the root list of the heap
  // If there is no m_child list, then do no work
  
  if (Result->m_child == NULL) {
    if (m_min_root == Result) m_min_root = NULL;
  }
  
  // If m_min_root==Result then there was only one root tree, so the
  // root list is simply the m_child list of that node (which is
  // NULL if this is the last node in the list)
  
  else if (m_min_root == Result) {
    m_min_root = Result->m_child;
  }
  
  // If m_min_root is different, then the m_child list is pushed into a
  // new temporary heap, which is then merged by Union() onto the
  // root list of this heap.
  
  else {
    m_childHeap = new FibHeap();
    m_childHeap->m_min_root = Result->m_child;
  }
  
  // Complete the disassociation of the Result node from the heap
  
  if (Result->m_child != NULL)
    Result->m_child->m_parent = NULL;
  Result->m_child = Result->m_parent = NULL;
  
  // If there was a m_child list, then we now merge it with the
  //	rest of the root list
  
  if (m_childHeap) makeUnion(m_childHeap);
  
  // Consolidate heap to find new minimum and do reorganize work
  
  if (m_min_root != NULL) consolidate();
  
  // Return the minimum node, which is now disassociated with the heap
  // It has m_left, m_right, m_parent, m_child, m_mark and m_degree cleared.
  
  return Result;
}

//===========================================================================
// DecreaseKey() - O(lg n) actual; O(1) amortized
//
// The O(lg n) actual cost stems from the fact that the depth, and
// therefore the number of ancestor m_parents, is bounded by O(lg n).
//===========================================================================

int  FibHeap::decreaseKey(FibHeapNode *theNode, FibHeapNode& NewKey)
{
  FibHeapNode *them_parent;
  
  if (theNode==NULL || *theNode < NewKey) return NOTOK;
  
  *theNode = NewKey;
  
  them_parent = theNode->m_parent;
  if (them_parent != NULL && *theNode < *them_parent) {
    cut(theNode, them_parent);
    cascadingCut(them_parent);
  }
  
  if (*theNode < *m_min_root) m_min_root = theNode;
  
  return OK;
}

//===========================================================================
// Delete() - O(lg n) amortized; extractMin() dominates
//
// Notice that if we don't own the heap nodes, then we clear the
// m_neg_infinity_flag on the deleted node.  Presumably, the programmer
// will be reusing the node.
//===========================================================================

int  FibHeap::remove(FibHeapNode *theNode)
{
  FibHeapNode Temp;
  int Result;
  
  if (theNode == NULL) return NOTOK;
  
  Temp.m_neg_infinity_flag = 1;
  Result = decreaseKey(theNode, Temp);
  
  if (Result == OK)
    if (extractMin() == NULL) Result = NOTOK;
    
    if (Result == OK) {
      if (getHeapOwnership()) delete theNode;
      else theNode->m_neg_infinity_flag = 0;
    }
    
    return Result;
}

//========================================================================
// Print()
//
// Used internally for debugging purposes.  The function prints the key
// value for each node along the root list, then it calls itself on each
// m_child list.   
//========================================================================

void FibHeap::print(FibHeapNode *Tree, FibHeapNode *them_parent)
{
  FibHeapNode* Temp = NULL;
  
  if (Tree == NULL) Tree = m_min_root;
  
  Temp = Tree;
  do {
    if (Temp->m_left == NULL) Rprintf( "(m_left is NULL)" );
    Temp->print();
    if (Temp->m_parent != them_parent)
      Rprintf("(m_parent is incorrect)" );
    if (Temp->m_right == NULL)
      Rprintf( "(m_right is NULL)" );
    else if (Temp->m_right->m_left != Temp)
      Rprintf( "(Error in m_left link m_left) ->" );
    else Rprintf( " <-> " );
    
    Temp = Temp->m_right;
    
    //		if (kbhit() && getch() == 27) {
    //			cout << "Hit a key to resume or ESC to break\n";
    //			if (getch() == 27) break;
    //		}
  } while (Temp != NULL && Temp != Tree);
  Rprintf( "\n" );
  
  Temp = Tree;
  do {
    Rprintf( "m_children of " );
    Temp->print();
    Rprintf( ": " );
    if (Temp->m_child == NULL) Rprintf( "NONE\n" );
    else print(Temp->m_child, Temp);
    Temp = Temp->m_right;
  } while (Temp!=NULL && Temp != Tree);
  
  if (them_parent == NULL) {
    //		char ch;
    
    Rprintf( "\n\n\n" );
    //		cin >> ch;
  }
}

//===========================================================================
//===========================================================================

void FibHeap::exchange(FibHeapNode*& N1, FibHeapNode*& N2)
{
  FibHeapNode *Temp;
  
  Temp = N1;
  N1 = N2;
  N2 = Temp;
}

//===========================================================================
// Consolidate()
//
// Internal function that reorganizes heap as part of an extractMin().
// We must find new minimum in heap, which could be anywhere along the
// root list.  The search could be O(n) since all nodes could be on
// the root list.  So, we reorganize the tree into more of a binomial forest
// structure, and then find the new minimum on the consolidated O(lg n) sized
// root list, and in the process set each m_parent pointer to NULL, and count
// the number of resulting subtrees.
//
// Note that after a list of n inserts, there will be n nodes on the root
// list, so the first extractMin() will be O(n) regardless of whether or
// not we consolidate.  However, the consolidation causes subsequent
// extractMin() operations to be O(lg n).  Furthermore, the extra cost of
// the first extractMin() is covered by the higher amortized cost of
// Insert(), which is the real governing factor in how costly the first
// extractMin() will be.
//===========================================================================

void FibHeap::consolidate()
{
  FibHeapNode *x, *y, *w;
  FibHeapNode *A[1+8*sizeof(long)]; // 1+lg(n)
  int  I=0, Dn = 1+8*sizeof(long);
  short d;
  
  // Initialize the consolidation detection array
  
  for (I=0; I < Dn; I++) A[I] = NULL;
  
  // We need to loop through all elements on root list.
  // When a collision of degree is found, the two trees
  // are consolidated in favor of the one with the lesser
  // element key value.  We first need to break the circle
  // so that we can have a stopping condition (we can't go
  // around until we reach the tree we started with
  // because all root trees are subject to becoming a
  // m_child during the consolidation).
  
  m_min_root->m_left->m_right = NULL;
  m_min_root->m_left = NULL;
  w = m_min_root;
  
  do {
    //cout << "Top of Consolidate's loop\n";
    //Print(w);
    
    x = w;
    d = x->m_degree;
    w = w->m_right;
    
    // We need another loop here because the consolidated result
    // may collide with another large tree on the root list.
    
    while (A[d] != NULL) {
      y = A[d];
      if (*y < *x) exchange(x, y);
      if (w == y) w = y->m_right;
      link(y, x);
      A[d] = NULL;
      d++;
      //cout << "After a round of Linking\n";
      //Print(x);
    }
    A[d] = x;
    
  } while (w != NULL);
  
  // Now we rebuild the root list, find the new minimum,
  // set all root list nodes' m_parent pointers to NULL and
  // count the number of subtrees.
  
  m_min_root = NULL;
  m_num_trees = 0;
  for (I = 0; I < Dn; I++)
    if (A[I] != NULL) addToRootList(A[I]);
}

//===========================================================================
// The node y is removed from the root list and becomes a subtree of node x.
//===========================================================================

void FibHeap::link(FibHeapNode *y, FibHeapNode *x)
{
  // Remove node y from root list
  
  if (y->m_right != NULL) y->m_right->m_left = y->m_left;
  if (y->m_left != NULL) y->m_left->m_right = y->m_right;
  m_num_trees--;
  
  // Make node y a singleton circular list with a m_parent of x
  
  y->m_left = y->m_right = y;
  y->m_parent = x;
  
  // If node x has no m_children, then list y is its new m_child list
  
  if (x->m_child == NULL) {
    x->m_child = y;
  }
  
  // Otherwise, node y must be added to node x's m_child list
  
  else {
    y->m_left = x->m_child;
    y->m_right = x->m_child->m_right;
    x->m_child->m_right = y;
    y->m_right->m_left = y;
  }
  
  // Increase the degree of node x because it's now a bigger tree
  
  x->m_degree ++;
  
  // Node y has just been made a m_child, so clear its mark
  
  if (y->m_mark) m_num_marked_nodes--;
  y->m_mark = 0;
}

//===========================================================================
//===========================================================================

void FibHeap::addToRootList(FibHeapNode *x)
{
  if (x->m_mark) m_num_marked_nodes --;
  x->m_mark = 0;
  
  m_num_nodes--;
  insert(x);
}

//===========================================================================
// Remove node x from the m_child list of its m_parent node y
//===========================================================================

void FibHeap::cut(FibHeapNode *x, FibHeapNode *y)
{
  if (y->m_child == x) y->m_child = x->m_right;
  if (y->m_child == x) y->m_child = NULL;
  
  y->m_degree --;
  
  x->m_left->m_right = x->m_right;
  x->m_right->m_left = x->m_left;
  
  addToRootList(x);
}

//===========================================================================
// Cuts each node in m_parent list, putting successive ancestor nodes on the
// root list until we either arrive at the root list or until we find an
// ancestor that is unmarked.  When a mark is set (which only happens during
// a cascading cut), it means that one m_child subtree has been lost; if a
// second subtree is lost later during another cascading cut, then we move
// the node to the root list so that it can be re-balanced on the next
// consolidate. 
//===========================================================================

void FibHeap::cascadingCut(FibHeapNode *y)
{
  FibHeapNode *z = y->m_parent;
  
  while (z != NULL) {
    if (y->m_mark == 0) {
      y->m_mark = 1;
      m_num_marked_nodes++;
      z = NULL;
    }
    else {
      cut(y, z);
      y = z;
      z = y->m_parent;
    }
  }
}
