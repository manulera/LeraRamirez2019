// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef LATTICE_H
#define LATTICE_H

#include <cmath>
#include <iostream>
#include "assert_macro.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "real.h"

class Fiber;
class Simul;

/// Array of discrete binding sites
/** 
 A Lattice can have different purposes:
 - to limit the local density of attached Hand,
 - to simulate traffic-jam of motors moving on Fibers,
 - to simulate binding/unbinding of a continuous substance.
 .
 
 The first 2 points are implemented in a specialized Hand: Digit.
 
 The index of a site can be positive or negative, and corresponds
 to the abscissa along the Fiber:
 @code
 int index(real abscissa) { return (int)(round(abscissa/laUnit)); }
 @endcode
 
 The Lattice hold a 'number of molecule' at each site. 
 The linear density can be derived by dividing by the unit length:
 @code
 density = value(index) / unit();
 @endcode
 */
template <typename SITE>
class Lattice
{
public:

    /// type for a Lattice cell
    typedef SITE    site_type;
    
private:
        
    /// lowest valid index (can be negative or positive)
    int           laInf;
    
    /// highest valid index plus one (laSup > laInf)
    int           laSup;
    
    /// Original allocated memory ( laSite = laSite0 - laInf )
    site_type *   laSite0;
    
    /// distance between adjacent sites
    real          laUnit;
    
    /// Array of sites, of valid range [laInf, laSup[
    site_type *   laSite;
    
    //------------------------------------------------------------------------------
    #pragma mark - Allocation
    
    /// allocate sites for indices within [inf, sup[, conserving existing indices
    int allocateLattice(int inf, int sup, int margin=8)
    {
        assert_true( inf <= sup );
        
        if ( laSite )
        {
            // check if existing boundaries are sufficient
            if ( laInf <= inf  &&  sup <= laSup )
                return 0;
        }
        
        // add security margin:
        inf -= margin;
        sup += margin;
        
        if ( laSite )
        {
            // only extend the current coverage of the lattice:
            if ( laInf < inf ) inf = laInf;
            if ( laSup > sup ) sup = laSup;
        }
        
        //std::clog << this << " Lattice::allocate [" << inf << ", " << sup << "[\n";
        
        const int dim = sup - inf;
        
        site_type * laSite0_new = new site_type[dim];
        site_type * laSite_new  = laSite0_new - inf;
        
        // reset new array:
        for ( int s = 0; s < dim; ++s )
            laSite0_new[s] = 0;
        
        // transfer the current information, from the overlapping zone
        if ( laSite )
        {
            for ( int s = laInf; s < laSup; ++s )
                laSite_new[s] = laSite[s];
            delete[] laSite0;
        }
        
        laInf = inf;
        laSup = sup;
        laSite0 = laSite0_new;
        laSite  = laSite_new;
        
        //std::clog<<"Lattice "<<this<<" covers ["<<laInf<<","<<laSup<<"[\n";
        return dim;
    }
    
    
    /// free all occupied memory
    void deallocate()
    {
        //printf("Lattice %p realeased\n", this);
        if ( laSite0 )
            delete[] laSite0;
        laSite0 = 0;
        laSite  = 0;
        laInf   = 0;
        laSup   = 0;
    }
    
    /// 'intercept' access occuring via automatic type conversion:
    template<typename T> site_type   value(const T s) const;
    template<typename T> site_type&  site(const T s);
    template<typename T> site_type&  operator[](const T s);
    
    /// disable assignment operator
    Lattice & operator =(const Lattice &);
    
    //------------------------------------------------------------------------------
    #pragma mark - Construction
    
public:
    
    /// Constructor
    Lattice(real unit=-1)
    {
        laUnit  = unit;
        laSite0 = 0;
        laSite  = 0;
        laInf   = 0;
        laSup   = 0;
    }
    
    
    /// Copy constructor
    Lattice(const Lattice & lat)
    {
        //std::clog << this << " Lattice::Lattice [" << laInf << ", " << laSup << "[\n";
        
        //copy member variables:
        memcpy(this, &lat, sizeof(Lattice));
        
        // allocate memory:
        laSite0 = new site_type[laSup-laInf];
        laSite  = laSite0 - laInf;
        
        // copy information:
        for ( int s = laInf; s < laSup; ++s )
            laSite[s] = lat.laSite[s];
    }
    
    
    /// Destructor
    ~Lattice()    { deallocate(); }
    
    
    /// change distance betwen adjacent sites
    void unit(real u)
    {
        if ( laUnit != u )
        {
            //std::clog << "changing Lattice:unit from " << laUnit << " to " << u << "\n";
            //preserve the same abscissa range, with the new lattice unit
            real i = abscissa(laInf);
            real s = abscissa(laSup+1.0);
            laUnit = u;
            if ( laSite )
                allocateLattice(index(i), index(s)+1);
        }
    }
    
    /// set the range of valid abscissa
    void setRange(real i, real s)
    {
        if ( laUnit <= 0 )
            throw InvalidParameter("the lattice:unit was not set");
        
        ///\todo use special Lattice value to indicate the ends of the fibers
        allocateLattice(index(i)-3, index(s)+4);
    }
    
    //------------------------------------------------------------------------------
    #pragma mark - Index / Abscissa
    
    /// true  if memory was allocated
    bool          ready() const { return laUnit > 0 && laSite0; }

    /// distance between adjacent sites
    real          unit() const { return laUnit; }
    
    /// index of the site containing abscissa `a`
    int           index(const real a) const { return (int)floor(a/laUnit); }
    
    /// true if index 'i' is covered by the lattice
    bool          valid(const int i) const { return ( laInf <= i  &&  i < laSup ); }
    
    /// true if index 'i' is not covered by the lattice
    bool          invalid(const int i) const { return ( i < laInf  || laSup <= i ); }
    
    
    /// the site of index `h` covers the abscissa range `unit * h < s < unit * ( h + 1 )`
    /**
     abscissa(h) returns the abscissa corresponding to the beginning of the site.
     The range covered by site 'h' is [ abscissa(h), abscissa(h+1) ]
     */
    real          abscissa(const real s) const { return s * laUnit; }
    
    /// true if abscissa `a` corresponds to a site that is covered by the lattice
    bool          within(const real a) const { return laInf * laUnit <= a  &&  a <= laSup * laUnit; }

    //------------------------------------------------------------------------------
    #pragma mark - Access

    /// value at index s
    site_type*    data()             const  { return laSite; }

    /// first valid index
    int           inf()              const  { return laInf; }
    
    /// last valid index plus one
    int           sup()              const  { return laSup; }

    /// value at index s
    site_type&    value(const int s) const  { assert_true(valid(s)); return laSite[s]; }
    
    /// reference to Site at index s
    site_type&    operator[](const int s)   { assert_true(valid(s)); return laSite[s]; }
    
    /// reference to Site at abscissa a : The argument is different than for operator []
    site_type&    site(const real a)        { int s=index(a); assert_true(valid(s)); return laSite[s]; }
   
#ifdef MANU_LATTICE
    /// ADDED BY MANU, sites with a value
    /// similar to inc, but at the right bit depth
    void        inc(const int s,const unsigned int val) const { assert_true(valid(s)); laSite[s]+= 1<<val;}
    /// similar to dec, but at the right bit depth
    void        dec(const int s,const unsigned int val) const { assert_true(valid(s)); laSite[s]-= 1<<val;}
    /// similar to vacant, but at the right bit depth
    bool        vacant(const int s, const unsigned int val){return !(laSite[s] & (1<<val));}
#else
    /// increment value at index s
    void          inc(const int s)          { assert_true(valid(s)); ++laSite[s];}
    
    /// decrement value at index s
    void          dec(const int s)          { assert_true(valid(s)); assert_true(laSite[s]>0); --laSite[s];}
    /// true if index s is unoccupied
    bool          vacant(const int s) const { assert_true(valid(s)); return ( laSite[s]==0 ); }
    
#endif
    
    /// set all sites to `value`
    void clear(site_type value = 0)
    {
        for ( int s=laInf; s<laSup; ++s ) 
            laSite[s] = value;
    }


    
    //------------------------------------------------------------------------------
    #pragma mark - Transfer
    
    /// transfer lat->sites within `[s, e[` to *this, and set lat->sites to zero
    void take(Lattice<SITE> * lat, int is, int ie)
    {
        // check all boundaries
        if ( is < laInf )      is = laInf;
        if ( is < lat->laInf ) is = lat->laInf;
        
        if ( ie > laSup )      ie = laSup;
        if ( ie > lat->laSup ) ie = lat->laSup;
        
        SITE * s = laSite + is;
        SITE * e = laSite + ie;        
        SITE * o = lat->laSite + is;

        // transfer values
        while ( s < e )
        {
            *s = *o;
            *o = 0;
            ++s;
            ++o;
        }
    }
    
    
    /// sum all sites in [inf, e[; set sites to zero and return total in `res`
    void takeM(Lattice<SITE> * lat, const int e)
    {
        take(lat, laInf, e);
    }
    
    /// sum all sites in [s, sup]; set sites to zero and return total in `res`
    void takeP(Lattice<SITE> * lat, const int s)
    {
        take(lat, s, laSup);
    }
    
    //------------------------------------------------------------------------------
    #pragma mark - Collect

    /// sum all sites in `[s, e[`, set them to zero and return total in `res`
    template <typename SUM>
    void collectR(SUM& res, int is, int ie)
    {
        res = 0;
#if ( 1 )
        // check boundaries
        if ( is < laInf )
            is = laInf;
        
        if ( ie > laSup )
            ie = laSup;
        
        // collect
        for ( ; is < ie; ++is )
        {
            res += laSite[is];
            laSite[is] = 0;
        }
#else
        // check boundaries
        SITE * s = laSite + ( ( is < laInf ) ? laInf : is );
        SITE * e = laSite + ( ( ie > laSup ) ? laSup : ie );        
        
        // collect
        while ( s < e )
        {
            res += *s;
            *s = 0;
            ++s;
        }
#endif
    }
    
    
    /// sum all sites in [inf, e[; set sites to zero and return total in `res`
    template <typename SUM>
    void collectM(SUM& res, const int e)
    {
        collectR(res, laInf, e);
    }
    
    /// sum all sites in ]s, sup]; set sites to zero and return total in `res`
    template <typename SUM>
    void collectP(SUM& res, const int s)
    {
        collectR(res, s+1, laSup);
    }
    
    /// sum of values in the entire lattice; set sites to zero
    template <typename SUM>
    void collect(SUM& res) const
    {
        collectR(res, laInf, laSup);
    }
    
    //------------------------------------------------------------------------------
    #pragma mark - Sums

    /// sum of values for all sites in `[s, e[` and return total in `res`
    template <typename SUM>
    void sum(SUM& res, int is, int ie) const
    {
        // check boundaries
        SITE * s = laSite + ( ( is < laInf ) ? laInf : is );
        SITE * e = laSite + ( ( ie > laSup ) ? laSup : ie );        
        
        // collect
        res = 0;       
        while ( s < e )
        {
            res += *s;
            ++s;
        }
    }
    
    /// sum of values for all sites below `e` (not-included)
    template <typename SUM>
    void sumM(SUM& res, const int e) const
    {
        sum(res, laInf, e);
    }

    /// sum of values for all sites above `s` (included)
    template <typename SUM>
    void sumP(SUM& res, const int s) const
    {
        sum(res, s, laSup);
    }

    
    /// sum of values in the entire lattice
    template <typename SUM>
    void sum(SUM& res) const
    {
        sum(res, laInf, laSup);
    }

    /// sum of values in the entire lattice using 'real' as accumulator
    real sum() const
    {
        real res;
        sum(res, laInf, laSup);
        return res;
    }

    //------------------------------------------------------------------------------
    #pragma mark - I/O
    
    /// write to file
    void write(Outputter& out) const
    {
        const size_t siz = sizeof(site_type);
        out.writeInt32(laInf);
        out.writeInt32(laSup);
        out.writeFloat(laUnit);
        out.writeUInt32(siz);
        
        if ( laSite )
        {
            if ( siz == 1 )
            {
                for ( int s = laInf; s < laSup; ++s )
                    out.writeUInt8(laSite[s]);
            }
            else if ( siz == 2 )
            {
                for ( int s = laInf; s < laSup; ++s )
                    out.writeUInt16(laSite[s]);
            }
            else
            {
                for ( int s = laInf; s < laSup; ++s )
                    out.writeFloat(laSite[s]);
            }
        }
    }
    
    
    /// read from file
    void read(Inputter& in)
    {
        try {
            
            int inf = in.readInt32();
            int sup = in.readInt32();
            real un = in.readFloat();
            int siz = in.readUInt32();
            
            if ( sup <= inf )
                throw InvalidIO("incoherent Lattice dimensions");
            if ( un <= 0 )
                throw InvalidIO("incoherent Lattice unit");

            // adjust unit size:
            laUnit = un;

            allocateLattice(inf, sup);
            clear();
            
            if ( siz == 1 )
            {
                for ( int s = inf; s < sup; ++s )
                    laSite[s] = in.readUInt8();
            }
            else if ( siz == 2 )
            {
                for ( int s = inf; s < sup; ++s )
                    laSite[s] = in.readUInt16();
            }
            else
            {
                for ( int s = inf; s < sup; ++s )
                    laSite[s] = in.readFloat();
            }
            
                
        }
        catch( Exception & e ) {
            e << ", in Lattice::read()";
            throw;
        }
    }
    
    
    
    /// debug function
    int bad()
    {
        if ( laSite ) {
            if ( laInf > laSup )
                return 1;
        }
        return 0;
    }
    
};

#endif
