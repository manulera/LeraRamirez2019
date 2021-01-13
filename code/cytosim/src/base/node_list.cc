// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// doubly linked list, STL style, with acces by iterators,
// some additions to manipulate the list: sorting, unsorting, etc.

#include "node_list.h"
#include "assert_macro.h"
#include "random.h"


void NodeList::push_front(Node * n)
{
    //MSG("NodeList: push_front   %p in   %p\n", n, this);

    n->nPrev = 0;
    n->nNext = nFirst;
    if ( nFirst )
        nFirst->nPrev = n;
    else
        nLast = n;
    nFirst = n;
    ++nSize;
    assert_false(bad());
}


void NodeList::push_back(Node * n)
{
    //MSG("NodeList: push_back   %p in   %p\n", n, this);
    
    n->nPrev = nLast;
    n->nNext = 0;
    if ( nLast )
        nLast->nNext = n;
    else
        nFirst = n;
    nLast = n;
    ++nSize;
    assert_false(bad());
}


/**
 Transfer objects in `list` to the end of `this`, until `list` is empty.
 */
void NodeList::merge(NodeList& list)
{
    Node * n = list.nFirst;
    
    if ( n )
    {
        if ( nLast )
            nLast->nNext = n;
        else
            nFirst = n;
        
        n->nPrev = nLast;
        nLast = list.nLast;
        nSize += list.nSize;
        
        list.nSize  = 0;
        list.nFirst = 0;
        list.nLast  = 0;
    }
}


void NodeList::push_after(Node * p, Node * n)
{
    n->nPrev = p;
    n->nNext = p->nNext;
    if ( p->nNext )
        p->nNext->nPrev = n;
    else
        nLast = n;
    p->nNext = n;
    ++nSize;
}


void NodeList::push_before(Node * p, Node * n)
{
    n->nNext = p;
    n->nPrev = p->nPrev;
    if ( p->nPrev )
        p->nPrev->nNext = n;
    else
        nFirst = n;
    p->nPrev = n;
    ++nSize;
}


void NodeList::pop_front()
{
    assert_true( nFirst );
 
    nFirst = nFirst->nNext;
    if ( nFirst )
        nFirst->nPrev = 0;
    else
        nLast = 0;
    --nSize;
}


void NodeList::pop_back()
{
    assert_true( nLast );
    
    nLast = nLast->nPrev;
    if ( nLast )
        nLast->nNext = 0;
    else
        nFirst = 0;
    --nSize;
}


void NodeList::pop(Node * n)
{
    assert_true( nSize > 0 );

    if ( n->nPrev )
        n->nPrev->nNext = n->nNext;
    else {
        assert_true( nFirst == n );
        nFirst = n->nNext;
    }
    
    if ( n->nNext )
        n->nNext->nPrev = n->nPrev;
    else {
        assert_true( nLast == n );
        nLast = n->nPrev;
    }
    
    n->nPrev = 0;
    n->nNext = 0;
    --nSize;
    assert_false(bad());
}


void NodeList::clear()
{
    Node * p, * n = nFirst;
    while ( n )
    {
        n->nPrev = 0;
        p = n->nNext;
        n->nNext = 0;
        n = p;
    }
    nFirst = 0;
    nLast  = 0;
    nSize  = 0;
}


void NodeList::erase()
{
    Node * n = nFirst;
    Node * p;
    while ( n )
    {
        p = n->nNext;
        delete(n);
        n = p;
    }
    nFirst = 0;
    nLast  = 0;
    nSize  = 0;
}


void NodeList::sort(int (*comp)(Node const*, Node const*))
{
    Node * ii = front();
    
    if ( ii == 0 )
        return;
    
    ii = ii->next();
        
    while ( ii )
    {
        Node * kk = ii->next();
        Node * jj = ii->prev();
            
        if ( comp(ii, jj) == 1 )
        {
            jj = jj->prev();
            
            while ( jj && comp(ii, jj) == 1 )
                jj = jj->prev();
            
            pop(ii);
            
            if ( jj )
                push_after(jj, ii);
            else
                push_front(ii);
        }
        ii = kk;
    }
}

/**
 Rearrange [F--P][Q--L] as [Q--L][F--P]
 */
void NodeList::permute(Node * p)
{
    if ( p  &&  p->nNext )
    {
        nLast->nNext   = nFirst;
        nFirst->nPrev  = nLast;
        nFirst         = p->nNext;
        nLast          = p;
        nLast->nNext   = 0;
        nFirst->nPrev  = 0;
    }
    assert_false( bad() );
}


/**
 Rearrange [F--P][X--Y][Q--L] as [X--Y][F--P][Q--L]
 
 If Q is between nFirst and P, this will destroy the list,
 but there is no way to check such condition here.
 */
void NodeList::shuffle_up(Node * p, Node * q)
{
    assert_true( p  &&  p->nNext );
    assert_true( q  &&  q->nPrev );
    
    if ( q != p->nNext )
    {
        nFirst->nPrev   = q->nPrev;
        q->nPrev->nNext = nFirst;
        nFirst          = p->nNext;
        nFirst->nPrev   = 0;
        p->nNext        = q;
        q->nPrev        = p;
    }
    assert_false( bad() );
}

/**
 Rearrange [F--P][X--Y][Q--L] as [F--P][Q--L][X--Y]
 
 If Q is between nFirst and P, this will destroy the list,
 but there is no way to check such condition here.
 */
void NodeList::shuffle_down(Node * p, Node * q)
{
    assert_true( p  &&  p->nNext );
    assert_true( q  &&  q->nPrev );
    
    if ( q != p->nNext )
    {
        nLast->nNext    = p->nNext;
        p->nNext->nPrev = nLast;
        p->nNext        = q;
        nLast           = q->nPrev;
        nLast->nNext    = 0;
        q->nPrev        = p;
    }
    assert_false( bad() );
}


void NodeList::mix(Random& rng)
{
    if ( nSize < 2 )
        return;
    
    unsigned int pp = rng.pint(nSize);
    unsigned int qq = rng.pint(nSize);

    unsigned int n = 0;
    Node *p = nFirst, *q;
    
    if ( pp+1 < qq )
    {
        for ( ; n < pp; ++n )
            p = p->nNext;
        for ( q = p; n < qq; ++n )
            q = q->nNext;
        
        shuffle_up(p, q);
    }
    else if ( qq+1 < pp )
    {
        for ( ; n < qq; ++n )
            p = p->nNext;
        for ( q = p; n < pp; ++n )
            q = q->nNext;
        
        shuffle_down(p, q);
    }
    else
    {
        for ( ; n < qq; ++n )
            p = p->nNext;

        permute(p);
    }
    assert_false(bad());
}


void NodeList::mix3(Random& rng)
{
    mix(rng);
    mix(rng);
    mix(rng);
}


unsigned int NodeList::count() const
{
    unsigned int cnt = 0;
    Node * p = nFirst;
    while ( p )
    {
        ++cnt;
        p = p->nNext;
    }
    return cnt;
}



bool NodeList::check(Node const* n) const
{
    Node * p = nFirst;
    while ( p )
    {
        if ( p == n )
            return true;
        p = p->nNext;
    }
    return false;
}



int NodeList::bad() const
{
    unsigned int cnt = 0;
    Node * p = nFirst, * q;
    
    if ( p  &&  p->nPrev != 0 )
        return 71;
    while ( p )
    {
        q = p->nNext;
        if ( q == 0 ) {
            if ( p != nLast )
                return 73;
        }
        else {
            if ( q->nPrev != p )
                return 74;
        }
        p = q;
        ++cnt;
    }
    
    if ( cnt != nSize )
        return 75;
    return 0;
}


