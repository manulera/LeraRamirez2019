// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef ARRAY_H
#define ARRAY_H

#include "assert_macro.h"
#include "random.h"

/** 
 Array<typename VAL> stores instatiations of the class VAL.
 VAL need to have a public constructor without argument.
 
 This class is similar to std::array<VAL> and resembles std::vector<VAL>.
 Many functions of std::vector are missing, 
 and some are original: remove_pack(), sort() and mix().
 
 Some functions:
 - allocate(int s) ensures that array[s-1] can then be accessed, as with C-arrays. 
 - operator[](int index) returns the object stored at that index.
 - mix() permutes the values to produce a random ordering.
 - data() allows direct access to the underlying C-array.
 .
 
 New memory is allocated if necessary by allocate(), and the values
 from the old array are copied to the new memory space.
 
 Allocation when it is done exceeds what is necessary by a bit,
 to ensure that allocation only occurs from time-to-time,
 even if one adds objects one by one to the array.
 
 Destruction of objects is not handled by the Array.
 */

/// Dynamic array of VAL
template <typename VAL>
class Array
{
public:

    /// typedef for the template argument
    typedef VAL value_type;

    /// iterator class type
    typedef value_type * iterator;
    
    /// const iterator class type
    typedef value_type const* const_iterator;

private:
    
    /// C-array holding the VAL
    VAL * val_;
    
    /// size of memory that was allocated for val_[]
    size_t alc_;
    
    /// number of objects currently present in the array
    size_t nbo_;
    
    /// size of the chunk used for memory allocation (a power of 2)
    size_t chk_;
    
#pragma mark -
private:
    
    /// the integer above s that is a multiple of chk_
    size_t chunked(size_t s)
    {
        return ( s + chk_ - 1 ) & ~( chk_ - 1 );
    }
    
    /// return smallest power of 2 that is greater or equal to `s`
    size_t next_power(size_t s)
    {
        if ( s & (s-1) )
        {
            do
                s &= s-1;
            while ( s & (s-1) );
            
            s <<= 1;
        }
        return s;
    }
    
    /// copy data
    inline void copy(VAL * dst, VAL const* src, const size_t cnt)
    {
        for ( size_t n = 0; n < cnt; ++n )
            dst[n] = src[n];
    }
    
#pragma mark -
public:
        
    /// Default creator without allocation
    Array() : val_(0), alc_(0), nbo_(0), chk_(16)
    {
    }

    /// set chunk size to `k` (do not allocate)
    Array(size_t k) : val_(0), alc_(0), nbo_(0)
    {
        if ( !k )
        {
            fprintf(stderr, "Array::chunk must not be null");
            exit(1);
        }

        chk_ = next_power(k);
    }
    
    /// Copy constructor
    Array(Array<VAL> const & o) : val_(0), alc_(0), nbo_(o.nbo_), chk_(o.chk_)
    {
        if ( o.alc_ )
        {
            //printf("Array %p copy constructor size %i\n", this, nbo_);
            allocate(o.alc_);
            copy(val_, o.val_, o.alc_);
        }
    }

    /// Destructor
    virtual ~Array()
    {
        deallocate();
    }
    
    /// Assignment operator
    Array& operator =(Array<VAL> const & o)
    {
        if ( o.nbo_ > alc_ )
        {
            //printf("Array %p allocated %i = from size %i\n", this, alc_, o.nbo_);
            deallocate();
            allocate(o.nbo_);
        }
        nbo_ = o.nbo_;
        copy(val_, o.val_, nbo_);
        return *this;
    }
    
    
#pragma mark -
    /// Number of objects
    unsigned size() const
    {
        return nbo_;
    }

    /// true if this Array holds no value
    bool  empty() const
    {
        return ( nbo_ == 0 );
    }
    
    /// Currently allocated size
    size_t capacity() const
    {
        return alc_;
    }
    
    /// Address of the underlying C-array
    VAL * data()
    {
        return val_;
    }
    
    /// Address of the underlying C-array
    VAL const * data() const
    {
        return val_;
    }
    
    /// pointer to first element
    iterator begin() const
    {
        return val_;
    }
    
    /// pointer to a position just past the last element
    iterator end() const
    {
        return val_+nbo_;
    }
    
    /// reference to Object at index ii (val_[ii])
    VAL & at(const size_t ii) const
    {
        assert_true( ii < nbo_ );
        return val_[ii];
    }
    
    /// reference to Object at index ii (val_[ii])
    VAL & operator[](const size_t ii) const
    {
        assert_true( ii < nbo_ );
        return val_[ii];
    }
    
    
#pragma mark -
    /// Allocate to hold `s` objects: valid indices are 0 <= indx < max
    void reallocate(const size_t alc_new)
    {
        VAL * val_new = new VAL[alc_new];
        if ( val_ )
        {
            if ( alc_ < alc_new )
                copy(val_new, val_, alc_);
            else
                copy(val_new, val_, alc_new);
            delete[] val_;
        }
        alc_ = alc_new;
        val_ = val_new;
    }
    
    /// Allocate to hold at least `s` objects: valid indices are 0 <= indx < max
    size_t allocate(const size_t s)
    {
        if ( s > alc_ )
        {
            reallocate(chunked(s));
            assert_true( alc_ >= s );
            return s;
        }
        return 0;
    }
    
    /// Allocate and set new values to `zero`
    size_t allocate_zero(const size_t size, VAL const& zero)
    {
        size_t res = allocate(size);
        if ( res )
        {
            //set the newly allocated memory to zero
            for ( size_t ii = nbo_; ii < alc_; ++ii )
                val_[ii] = zero;
        }
        return res;
    }
    
    /// truncate Array to a smaller size
    void truncate(const size_t size)
    {
        if ( size < nbo_ )
            nbo_ = size;
    }
    
    /// Set the size of this Array to `size` (allocate or truncate if necessary)
    void resize(const size_t size)
    {
        if ( size < nbo_ )
            nbo_ = size;
        else if ( size > nbo_ )
        {
            allocate(size);
            nbo_ = size;
        }
    }
    
    /// Release occupied memory
    void deallocate()
    {
        if ( val_ )
        {
            //printf("Array %p deallocate %i\n", this, allocated);
            delete[] val_;
            val_ = 0;
        }
        alc_ = 0;
        nbo_ = 0;
    }
    
    /// Set the number of objects to zero
    inline void clear()
    {
        nbo_ = 0;
    }
    
    /// Delete all values as if they were pointers to Object
    void destroy()
    {
        assert_true( val_ || nbo_==0 );
        for ( size_t ii=0; ii < nbo_; ++ii )
        {
            if ( val_[ii] )
            {
                //std::clog << " delete " << val_[ii] << std::endl;
                delete(val_[ii]);
                val_[ii] = 0;
            }
        }
        nbo_ = 0;
    }
    
    /// Set all values to `zero`
    void zero(VAL const& zero)
    {
        assert_true( val_ || alc_==0 );
        for ( size_t ii=0; ii < alc_; ++ii )
            val_[ii] = zero;
    }
    
    
#pragma mark -
    
    /// Increment the size of the array, and return new value at end of it
    VAL & new_val()
    {
        if ( nbo_ >= alc_ )
            reallocate(chunked(nbo_+1));
        VAL& res = val_[nbo_++];
        return res;
    }
    
    /// Add `np` at the end of this Array
    void push_back(const VAL np)
    {
        if ( nbo_ >= alc_ )
            reallocate(chunked(nbo_+1));
        val_[nbo_++] = np;
    }
    
    /// Add the elements of `array` at the end of this Array
    void append(const Array<VAL> array)
    {
        allocate(nbo_+array.nbo_);
        for ( size_t ii = 0; ii < array.nbo_; ++ii )
            val_[ii+nbo_] = array.val_[ii];
        nbo_ += array.nbo_;
    }
    
    /// Add the elements of `array` at the end of this Array
    void append_except(const Array<VAL> array, const VAL& val)
    {
        allocate(nbo_+array.nbo_);
        for ( size_t ii = 0; ii < array.nbo_; ++ii )
            if ( array.val_[ii] != val )
                val_[nbo_++] = array.val_[ii];
    }

    /// Return index of `obj`, or -1 if not found in the list (linear search)
    int index(const VAL obj) const
    {
        assert_true( val_ || nbo_==0 );
        for ( size_t ii = 0; ii < nbo_; ++ii )
            if ( val_[ii] == obj )
                return ii;
        return -1;
    }
    
    /// Replace `old_value` by `new_value`, or return false if `old_value` is not found
    bool replace(VAL const& old_value, VAL const& new_value)
    {
        assert_true( val_ || nbo_==0 );
        for ( size_t ii=0; ii < nbo_; ++ii )
        {
            if ( val_[ii] == old_value )
            {
                val_[ii] = new_value;
                return true;
            }
        }
        return false;
    }
    
#pragma mark -
    
    /// set `np` at any position equal to `zero`, or at the end of the array
    void push_pack(const VAL np, const VAL& zero)
    {
        assert_true( val_ || nbo_==0 );
        for ( size_t ii = 0; ii < nbo_; ++ii )
        {
            if ( val_[ii] == zero )
            {
                val_[ii] = np;
                return;
            }
        }
        push_back(np);
    }

    
    /// Returns the number of occurence of 'val' in the array
    size_t count(VAL const& val) const
    {
        if ( val_ == 0 || nbo_==0 )
            return 0;
        size_t res = 0;
        for ( size_t ii = 0; ii < nbo_; ++ii )
            if ( val_[ii] == val ) ++res;
        return res;
    }

    /// Number of values which are different from `val` in the array
    size_t count_except(VAL const& val) const
    {
        if ( val_ == 0 || nbo_==0 )
            return 0;
        size_t res = 0;
        for ( size_t ii = 0; ii < nbo_; ++ii )
            if ( val_[ii] != val ) ++res;
        return res;
    }

    
    /**
     Remove all entries which are equal to `zero`, and pack array by shuffling values around.
     The order of the elements is not preserved, and copy operations are minimized
     */
    template <typename T>
    static T * remove_pack(T * s, T * e, T const& zero)
    {
        if ( e <= s )
            return s;
        --e;
        while ( s < e )
        {
            // find the next `zero` going upward:
            while ( *s != zero )
            {
                ++s;
                if ( e <= s )
                    return e + ( *e != zero );
            }
            // going downward, skip `zero` values:
            while ( *e == zero )
            {
                --e;
                if ( e <= s )
                    return e;
            }
            // flip the two values:
            *s = *e;
            *e = zero; // maybe not necessary
            ++s;
            --e;
        }
        return e + ( *e != zero );
    }
    
    
    

    /// Remove all entries which are equal to `zero`, and pack array
    void remove_pack(VAL const& zero)
    {
        assert_true( val_ || nbo_==0 );
        nbo_ = remove_pack(val_, val_+nbo_, zero) - val_;
    }
    
    
    /// Sort array using `std::qsort()`
    void sort(int (*comp)(const void *, const void *))
    {
        if (val_)
            qsort(val_, nbo_, sizeof(VAL), comp);
    }
    
//    /// Custom sort created by manu
//    void sort(int (*comp)(const VAL *, const VAL *))
//    {
//        if (val_)
//            qsort(val_, nbo_, sizeof(VAL), comp);
//    }
    
    
    /// Return one of the value in the array, chosen randomly
    VAL& pick_one(Random& rng)
    {
        assert_true(nbo_>0);
        return val_[rng.pint(nbo_)];
    }

    /// Move the last Object on top, push all other values down by one slot
    void turn()
    {
        if ( nbo_ > 1 )
        {
            assert_true(val_);
        
            VAL * tmp = val_[0];
            for ( size_t ii = 0; ii < nbo_-1; ++ii )
                val_[ii] = val_[ii+1];
            val_[nbo_-1] = tmp;
        }
    }
    
    /// exchange the values of `a` and `b`
    static void swap(VAL& a, VAL& b)
    {
        VAL tmp = a;
        a = b;
        b = tmp;
    }
    
    /// Swap two random values in the array
    void permute(Random& rng)
    {
        assert_true(val_);
        size_t ii = rng.pint() % nbo_;
        size_t jj = rng.pint() % nbo_;
        if ( ii != jj )
            swap(val_[ii], val_[jj]);
    }
    
    
    /// Randomly permutes all objects in the array
    /**
     This produces uniform shuffling in linear time.
     see Knuth's The Art of Programming, Vol 2 chp. 3.4.2 
     */
    void mix(Random& rng)
    {
        assert_true( val_ || nbo_==0 );
        size_t jj = nbo_, kk;
        while ( jj > 1 )
        {
            kk = rng.pint() % jj;  //between 0 and j-1
            --jj;
            swap(val_[jj], val_[kk]);
        }
    }
};




#endif
