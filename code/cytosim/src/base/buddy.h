// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef BUDDY_H
#define BUDDY_H

#include <vector>
#include <algorithm>


/// Maintains mutual relationship between objects.
/**
 Buddy implements mutual relationship between objects.
 
 The class keeps track of a list of `buddies`.
 Relationship is established with hello().
 When an object is destroyed, it calls goodbye() for all its buddies.
 
 This class can be used when an object needs to know if another object is destroyed,
 and vice-versa.
 
 F. Nedelec 11 Aug. 2012
 */
class Buddy
{

private:
    
    /// type for a list of buddies
    typedef std::vector<Buddy *> BuddyList;
    
    /// list of buddies
    BuddyList buddies_;
    
private:
    
    /// add `b` into the list of buddies, or complain if already present
    void enroll(Buddy * b)
    {
#if ( 1 )
        // complain if buddy is known already:
        BuddyList::iterator bi = std::find(buddies_.begin(), buddies_.end(), b);
        if ( bi != buddies_.end() )
        {
            std::cerr << " Error: duplicated Buddy::enroll()" << std::endl;
            return;
        }
#endif
        
        // find an empty spot:
        bi = std::find(buddies_.begin(), buddies_.end(), (Buddy*)0);
        if ( bi != buddies_.end() )
            *bi = b;
        else
            buddies_.push_back(b);
    }
    
    
    /// replace the buddy that may have been at index `ix`
    void enroll(Buddy * b, unsigned ix)
    {
        if ( ix < buddies_.size() )
        {
            if ( buddies_[ix]  &&  buddies_[ix] != b )
            {
                goodbye(buddies_[ix]);
                buddies_[ix]->goodbye(this);
                buddies_[ix]->unlist(this);
            }
        }
        else
            buddies_.resize(ix+1, 0);
        
        buddies_[ix] = b;
    }

    
    /// removes `b` from the list of known buddy, do not call goodbye()
    Buddy * unlist(Buddy * b)
    {
        BuddyList::iterator bi = std::find(buddies_.begin(), buddies_.end(), b);
        if ( bi != buddies_.end() )
        {
            *bi = 0;
            return b;
        }
        return 0;
    }

public:
    
    /// constructor
    Buddy() {}
    
    /// upon destruction, goodbye is called for all buddies
    virtual ~Buddy()
    {
        for ( BuddyList::iterator bi = buddies_.begin(); bi < buddies_.end(); ++bi )
            if ( *bi )
            {
                goodbye(*bi);
                (*bi)->goodbye(this);
                (*bi)->unlist(this);
                *bi = 0;
            }
    }

    /// will make `this` and `guy` mutual buddies
    void connect(Buddy * guy)
    {
        if ( guy )
        {
            enroll(guy);
            guy->enroll(this);
        }
    }
    
    /// will the association, without calling goodbye()
    void disconnect(Buddy * guy)
    {
        if ( guy )
        {
            unlist(guy);
            guy->unlist(this);
        }
    }

    /// this is called everytime a known buddy is destroyed
    virtual void goodbye(Buddy *)
    {
    }
    

    /// returns the number of registered buddies
    unsigned nbBuddies() const
    {
        return buddies_.size();
    }
    
    
    /// return buddy at index `ix`
    Buddy * buddy(unsigned ix) const
    {
        if ( ix < buddies_.size() )
            return buddies_[ix];
        return 0;
    }
    
    
    /// returns true if `guy` is a buddy
    int check(Buddy * guy) const
    {
        BuddyList::const_iterator bi = std::find(buddies_.begin(), buddies_.end(), guy);
        
        return ( bi != buddies_.end() );
    }
    
};

#endif

