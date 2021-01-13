// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "vector3.h"


/**
 This accepts 'X Y Z' but also 'X' and 'X Y'.
 At least one scalar must be read to be valid
 */
std::istream & operator >> (std::istream& is, Vector3& v)
{
    if ( is >> v.XX )
    {
        if ( is >> v.YY )
        {
            if ( is >> v.ZZ )
                ;
            else
            {
                v.ZZ = 0;
                is.clear();
            }
        }
        else
        {
            v.YY = 0;
            v.ZZ = 0;
            is.clear();
        }
    }
    return is;
}


std::ostream & operator << (std::ostream& os, Vector3 const& v)
{
    std::streamsize w = os.width();
    os << v.XX << " ";
    os.width(w);
    os << v.YY << " ";
    os.width(w);
    os << v.ZZ;
    return os;
}


const Vector3 interpolate(const Vector3& a, real x, const Vector3& b)
{
    return Vector3(a.XX+x*b.XX, a.YY+x*b.YY, a.ZZ+x*b.ZZ);
}


real distanceSqr(const Vector3& a, const Vector3& b)
{
    return (a.XX-b.XX)*(a.XX-b.XX) + (a.YY-b.YY)*(a.YY-b.YY) + (a.ZZ-b.ZZ)*(a.ZZ-b.ZZ);
}


real distance(const Vector3& a, const Vector3& b)
{
    return sqrt((a.XX-b.XX)*(a.XX-b.XX) + (a.YY-b.YY)*(a.YY-b.YY) + (a.ZZ-b.ZZ)*(a.ZZ-b.ZZ));
}
