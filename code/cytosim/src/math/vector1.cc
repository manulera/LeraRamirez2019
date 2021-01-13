// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "vector1.h"

/**
 Optinally read two scalar values equal to zero
 */
void eatTwoZeros(std::istream& is)
{
    if ( is.good() )
    {
        std::streampos isp = is.tellg();
        real Z = 0;
        is >> Z;
        if ( is.good() && Z==0 )
        {
            isp = is.tellg();
            is >> Z;
            if ( is.good() && Z==0 )
                return;
        }
        // restore initial state:
        is.seekg(isp);
        is.clear();
    }
}


/**
 This accepts 'X 0 0' but also 'X' and 'X 0'.
 At least one scalar must be read to be valid
 */
std::istream & operator >> (std::istream& is, Vector1& v)
{
    if ( is >> v.XX )
        eatTwoZeros(is);
    return is;
}


std::ostream & operator << (std::ostream& os, Vector1 const& v)
{
    os << v.XX;
    return os;
}


const Vector1 interpolate(const Vector1& a, real x, const Vector1& b)
{
    return Vector1(a.XX+x*b.XX, 0);
}


real distanceSqr(const Vector1& a, const Vector1& b)
{
    return (a.XX-b.XX)*(a.XX-b.XX);
}


real distance(const Vector1& a, const Vector1& b)
{
    return fabs(a.XX-b.XX);
}

