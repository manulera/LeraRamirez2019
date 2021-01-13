// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef NODE_LIST_H
#define NODE_LIST_H

#include "node.h"
class Random;

/// Doubly linked list of Nodes
/**
 
 This class is similar to the standard template library <std::list>
 and the naming of the functions is consistent with STL whenever possible.
 
 The NodeList holds pointers to the first and last elements of the list,
 and it keeps track of the number of objects linked.
 Functions are given to link and unlink Nodes in constant time.\n
 
 A function mix() randomize the order of the Nodes in the list, which is
 necessary in a simulation to avoid any bias which could derive from fixed ordering.
 
 The list is zero-terminated on both sides, and it can be traversed in either ways:
 for ( Node * n = front(); n ; n = n->next() );
 for ( Node * n = back() ; n ; n = n->prev() );
 
 */

class NodeList
{
    
private:
        
    /// First Node of the list
    Node *          nFirst;
    
    /// Last Node of the list
    Node *          nLast;
    
    /// Number of Node in the list
    unsigned int    nSize;
    
    /// Disabled copy constructor
    NodeList(NodeList const&);
    
    /// Disabled copy assignment
    NodeList& operator=(NodeList const&);
    
public:
    
    /// default constructor
    NodeList() : nFirst(0), nLast(0), nSize(0) { }
    
    /// Destructor
    virtual         ~NodeList()         { clear(); }
    
    /// First Node in list
    Node *          front()       const { return nFirst; }
    
    /// Last Node in list
    Node *          back()        const { return nLast; }

    /// Number of objects in the list
    unsigned int    size()        const { return nSize; }
    
    /// true if list has zero elements
    bool            empty()       const { return nFirst == 0; }
    
    /// put Node first in the list
    void            push_front(Node *);
    
    /// put Node last in the list
    void            push_back(Node *);
    
    /// import all objects from given list, and empty it
    void            merge(NodeList& list);
    
    /// put new Node `n` after existing one `p`
    void            push_after(Node * p, Node * n);
    
    /// put new Node `n` before existing one `p`
    void            push_before(Node * p, Node * n);
    
    /// Remove Node `n` from list
    void            pop(Node * n);
    
    /// Remove top Node from list
    void            pop_front();

    /// Remove top Node from list
    void            pop_back();
   
    /// clear the list
    void            clear();
    
    /// delete all nodes, clearing the list on the way
    void            erase();
    
    /// sort according to given function
    void            sort(int (*comp)(Node const*, Node const*));
    
    /// Rearrange the list by exchanging the portions before and after `p`
    void            permute(Node * p);
    
    /// Rearrange the list by moving a central portion to the top
    void            shuffle_up(Node * p, Node * q);
    
    /// Rearrange the list by moving a central portion to the bottom
    void            shuffle_down(Node * p, Node * q);
    
    /// Mix list using permute() and shuffle() functions
    void            mix(Random&);
    
    /// call mix() three times
    void            mix3(Random&);

    /// count number of elements in the list
    unsigned int    count() const;
    
    /// return `true` if element appears in the list
    bool            check(Node const* n) const;
    
    /// test coherence of list
    int             bad() const;
};


#endif
